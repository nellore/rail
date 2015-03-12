#!/usr/bin/env bash
# Patches MultipartUploadManager in EMRFS so it tries to upload each part DEFAULT_PART_ATTEMPTS = 100 before failing
set -e

mkdir -p MultipartWork
cd MultipartWork
cat >MultipartUploadManager.java <<EOF
package com.amazon.ws.emr.hadoop.fs.s3;

import com.amazon.ws.emr.hadoop.fs.cse.CSEUtils;
import com.amazonaws.AmazonClientException;
import com.amazonaws.AbortedException;
import com.amazonaws.AmazonWebServiceRequest;
import com.amazonaws.event.ProgressListener;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.model.AbortMultipartUploadRequest;
import com.amazonaws.services.s3.model.CompleteMultipartUploadRequest;
import com.amazonaws.services.s3.model.CompleteMultipartUploadResult;
import com.amazonaws.services.s3.model.InitiateMultipartUploadRequest;
import com.amazonaws.services.s3.model.InitiateMultipartUploadResult;
import com.amazonaws.services.s3.model.ObjectMetadata;
import com.amazonaws.services.s3.model.PartETag;
import com.amazonaws.services.s3.model.UploadPartRequest;
import com.amazonaws.services.s3.model.UploadPartResult;
import com.google.common.base.Preconditions;
import com.google.common.base.Strings;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.util.concurrent.FutureCallback;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.Future;
import org.apache.hadoop.conf.Configuration;
import org.joda.time.DateTime;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MultipartUploadManager {
    public static final Logger LOG = LoggerFactory.getLogger((Class)MultipartUploadManager.class);
    private static final double DEFAULT_TH_FRACTION_PARTS_COMPLETED = 0.5;
    private static final double DEFAULT_FRACTION_PART_AVG_COMPLETION_TIME = 1.0;
    private static final int MIN_PART_ATTEMPTS = 2;
    private static final int DEFAULT_PART_ATTEMPTS = 5;
    private static final double TH_FRACTION_MAX_PART_SIZE = 0.7;
    private String bucketName;
    private String key;
    private int partNumber = 0;
    private ConcurrentLinkedQueue<Future<UploadPartResult>> parts = new ConcurrentLinkedQueue();
    private double thFractionPartsCompleted;
    private double fractionPartAvgCompletionTime;
    private int partAttempts;
    private Set<Integer> incompletePartNums = Collections.synchronizedSet(new HashSet());
    private ConcurrentHashMap<Integer, List<MultiPartUploadFuture>> partNumFutureMap = new ConcurrentHashMap();
    private Object partNumFutureMapHandle = new Object();
    private String uploadId;
    private ProgressListener progressListener;
    private Map<String, String> userMetadata = Maps.newHashMap();
    private Object desiredFilename = null;
    private String serverSideEncryptionAlgorithm;
    private AmazonS3 s3;
    private ListeningExecutorService executorService;
    private boolean closed = false;
    private long totalLength = 0;
    private long partSize;
    private Configuration conf;

    private MultipartUploadManager() {
    }

    private void initializeFromConf(Configuration conf) {
        this.thFractionPartsCompleted = conf.getDouble("fs.s3.multipart.th.fraction.parts.completed", 0.5);
        if (this.thFractionPartsCompleted <= 0.0 || this.thFractionPartsCompleted >= 1.0) {
            LOG.warn("The value of fs.s3.multipart.th.fraction.parts.completed is not in valid range: (0,1). Will set to default: 0.5");
            this.thFractionPartsCompleted = 0.5;
        }
        this.fractionPartAvgCompletionTime = conf.getDouble("fs.s3.multipart.fraction.part.avg.completion.time", 1.0);
        if (this.fractionPartAvgCompletionTime <= 0.0 || this.fractionPartAvgCompletionTime > 1.0) {
            LOG.warn("The value of fs.s3.multipart.fraction.part.avg.completion.time is not in valid range: (0,1]. Will set to default: 1.0");
            this.fractionPartAvgCompletionTime = 1.0;
        }
        this.partAttempts = conf.getInt("fs.s3.multipart.part.attempts", 5);
        if (this.partAttempts < 2) {
            LOG.warn("The value of fs.s3.multipart.part.attempts is less than min value: 2. Will default to min value.");
            this.partAttempts = 2;
        }
    }

    public double getThFractionPartsCompleted() {
        return this.thFractionPartsCompleted;
    }

    public double getFractionPartAvgCompletionTime() {
        return this.fractionPartAvgCompletionTime;
    }

    public int getPartAttempts() {
        return this.partAttempts;
    }

    public void start() {
        Preconditions.checkState((boolean)(this.partNumber == 0));
        this.partNumber = 1;
        InitiateMultipartUploadRequest request = new InitiateMultipartUploadRequest(this.bucketName, this.key);
        ObjectMetadata md = new ObjectMetadata();
        if (!Strings.isNullOrEmpty((String)this.serverSideEncryptionAlgorithm)) {
            md.setSSEAlgorithm(this.serverSideEncryptionAlgorithm);
        }
        md.setUserMetadata(this.userMetadata);
        md.setContentType("binary/octet-stream");
        if (this.desiredFilename != null) {
            md.setContentDisposition("filename=\"" + this.desiredFilename + "\"");
        }
        request.setObjectMetadata(md);
        InitiateMultipartUploadResult res = this.s3.initiateMultipartUpload(request);
        this.uploadId = res.getUploadId();
    }

    public void abort() throws IOException {
        AbortMultipartUploadRequest request = new AbortMultipartUploadRequest(this.bucketName, this.key, this.uploadId);
        try {
            for (Map.Entry<Integer, List<MultiPartUploadFuture>> entry : this.partNumFutureMap.entrySet()) {
                List<MultiPartUploadFuture> multiPartUploadFutures = entry.getValue();
                for (MultiPartUploadFuture multipartUploadFuture : multiPartUploadFutures) {
                    Future<UploadPartResult> future = multipartUploadFuture.getFuture();
                    future.cancel(true);
                }
            }
            this.s3.abortMultipartUpload(request);
            this.partNumber = 0;
            this.parts.clear();
            this.uploadId = null;
        }
        catch (AmazonClientException e) {
            throw new IOException((Throwable)e);
        }
    }

    private MultiPartUploadFuture createMultiPartUploadFuture(long partSize, MultipartUploadCallable multipartUploadCallable, int partNum) {
        ListenableFuture future = this.executorService.submit((Callable)multipartUploadCallable);
        MultipartUploadFutureCallBack futureCallBack = new MultipartUploadFutureCallBack(future, partNum);
        Futures.addCallback((ListenableFuture)future, (FutureCallback)futureCallBack);
        MultiPartUploadFuture multiPartUploadFuture = new MultiPartUploadFuture(partSize, (Future<UploadPartResult>)future, multipartUploadCallable);
        DateTime startTime = DateTime.now();
        multiPartUploadFuture.setStartTime(startTime);
        return multiPartUploadFuture;
    }

    private boolean shouldSpawnNewFuture(MultiPartUploadFuture remFuture) {
        long totalTime = 0;
        int completedFutures = 0;
        for (Map.Entry<Integer, List<MultiPartUploadFuture>> entry : this.partNumFutureMap.entrySet()) {
            MultiPartUploadFuture multipartUploadFuture;
            int partNum = entry.getKey();
            if (this.incompletePartNums.contains(partNum) || (double)(multipartUploadFuture = entry.getValue().get(0)).getPartSize() / (double)this.partSize < 0.7) continue;
            totalTime+=multipartUploadFuture.getEndTime().getMillis() - multipartUploadFuture.getStartTime().getMillis();
            ++completedFutures;
        }
        if (completedFutures == 0) {
            return false;
        }
        double avgCompletionTime = (double)totalTime / (double)completedFutures;
        long timeSinceFutureStarted = DateTime.now().getMillis() - remFuture.getStartTime().getMillis();
        double fractionPartsIncomplete = (double)this.incompletePartNums.size() / (double)this.partNumFutureMap.size();
        double thresholdTimeForSpawn = (1.0 + fractionPartsIncomplete) * (avgCompletionTime * this.fractionPartAvgCompletionTime);
        LOG.debug("Threshold time before spawn: " + thresholdTimeForSpawn);
        LOG.debug("Time since incomplete future started: " + timeSinceFutureStarted);
        if ((double)timeSinceFutureStarted > thresholdTimeForSpawn) {
            LOG.debug("Incomplete future exceeded threshold, will start new one..");
            return true;
        }
        return false;
    }

    private void spawnNewFutureIfNeeded(int partNum) throws IOException {
        List<MultiPartUploadFuture> multiPartUploadFutures = this.partNumFutureMap.get(partNum);
        LOG.debug("Number of running attempts for: " + partNum + " are: " + multiPartUploadFutures.size());
        MultiPartUploadFuture multiPartUploadFuture = multiPartUploadFutures.get(multiPartUploadFutures.size() - 1);
        if (this.shouldSpawnNewFuture(multiPartUploadFuture)) {
            if (multiPartUploadFutures.size() >= this.partAttempts) {
                LOG.error("Upload attempts for part num: " + partNum + " have already reached max limit of: " + this.partAttempts + ", will throw exception and fail");
                throw new IllegalStateException("Reached max limit of upload attempts for part");
            }
            LOG.debug("Creating new future for partNum: " + partNum);
            MultipartUploadCallable multipartUploadCallable = multiPartUploadFuture.getMultiPartUploadCallable();
            File origPartFile = multipartUploadCallable.getPartFile();
            File clonePartFile = File.createTempFile("emrfs-", ".tmp", origPartFile.getParentFile());
            try {
                Files.copy((File)origPartFile, (File)clonePartFile);
            }
            catch (IOException ioe) {
                if (!origPartFile.exists()) {
                    LOG.warn("Caught IOException while trying to create clone temp part file for partNum: " + partNum + ". Seems the other thread finished first and deleted the original temp part file", (Throwable)ioe);
                    return;
                }
                LOG.error("Caught IOException while trying to create clone temp part file for partNum: " + partNum + " " + ioe);
                throw ioe;
            }
            MultipartUploadCallable newMultiPartUploadCallable = new MultipartUploadCallable(multipartUploadCallable.getPartNumber(), clonePartFile, multipartUploadCallable.getStart());
            newMultiPartUploadCallable.setIsAnOriginalPartFile(false);
            MultiPartUploadFuture newMultipartUploadFuture = this.createMultiPartUploadFuture(multiPartUploadFuture.getPartSize(), newMultiPartUploadCallable, partNum);
            multiPartUploadFutures.add(newMultipartUploadFuture);
        }
    }

    public CompleteMultipartUploadResult commit() throws IOException {
        try {
            int totalParts = this.partNumFutureMap.size();
            while (this.incompletePartNums.size() != 0) {
                int completedFutures = totalParts - this.incompletePartNums.size();
                if (completedFutures > 0 && (double)completedFutures >= this.thFractionPartsCompleted * (double)totalParts) {
                    LOG.debug("" + completedFutures + " part(s) completed, checking heuristic...");
                    Set<Integer> set = this.incompletePartNums;
                    synchronized (set) {
                        Iterator<Integer> i$ = this.incompletePartNums.iterator();
                        while (i$.hasNext()) {
                            int partNum = i$.next();
                            this.spawnNewFutureIfNeeded(partNum);
                        }
                    }
                }
                Thread.sleep(1000);
            }
            ArrayList<PartETag> etags = new ArrayList<PartETag>();
            for (Map.Entry<Integer, List<MultiPartUploadFuture>> entry : this.partNumFutureMap.entrySet()) {
                MultiPartUploadFuture multiPartUploadFuture = entry.getValue().get(0);
                etags.add(multiPartUploadFuture.getFuture().get().getPartETag());
            }
            LOG.debug("Complete multipart upload " + this.uploadId + " with bucket '" + this.bucketName + "' key '" + this.key + "' and etags '" + etags + "'");
            CompleteMultipartUploadResult completeMultipartResult = this.s3.completeMultipartUpload(new CompleteMultipartUploadRequest(this.bucketName, this.key, this.uploadId, etags));
            LOG.info("completed multipart upload of " + this.partNumFutureMap.size() + " parts " + this.totalLength + " bytes");
            CSEUtils.deletePreviousInstructionFileIfNecessary((Configuration)this.conf, (AmazonS3)this.s3, (String)this.bucketName, (String)this.key);
            this.closed = true;
            return completeMultipartResult;
        }
        catch (Exception e) {
            LOG.info("completeMultipartUpload error for key: " + this.key, (Throwable)e);
            this.abort();
            throw new IOException("Error closing multipart upload", e);
        }
    }

    public void addPartFromFile(File contents, long start, long length) throws IOException {
        MultipartUploadCallable multipartUploadCallable = new MultipartUploadCallable(this.partNumber, contents, start);
        multipartUploadCallable.setIsAnOriginalPartFile(true);
        MultiPartUploadFuture multiPartUploadFuture = this.createMultiPartUploadFuture(length, multipartUploadCallable, this.partNumber);
        ArrayList<MultiPartUploadFuture> multiPartUploadFutures = new ArrayList<MultiPartUploadFuture>();
        multiPartUploadFutures.add(multiPartUploadFuture);
        Object object = this.partNumFutureMapHandle;
        synchronized (object) {
            this.partNumFutureMap.put(this.partNumber, multiPartUploadFutures);
            this.partNumFutureMapHandle.notify();
        }
        this.incompletePartNums.add(this.partNumber);
        this.totalLength+=contents.length();
        ++this.partNumber;
    }

    public int getPartNumber() {
        return this.partNumber;
    }

    public ConcurrentLinkedQueue<Future<UploadPartResult>> getParts() {
        return this.parts;
    }

    public int countPendingUploads() {
        int count = 0;
        for (Map.Entry<Integer, List<MultiPartUploadFuture>> entry : this.partNumFutureMap.entrySet()) {
            List<MultiPartUploadFuture> multiPartUploadFutures = entry.getValue();
            for (MultiPartUploadFuture multipartUploadFuture : multiPartUploadFutures) {
                Future<UploadPartResult> future = multipartUploadFuture.getFuture();
                if (future.isDone()) continue;
                ++count;
            }
        }
        return count;
    }

    private class MultipartUploadCallable
    implements Callable<UploadPartResult> {
        private final int partNumber;
        private final File partFile;
        private final long start;
        private boolean shouldCallAbortOnCompletion;
        private boolean isAnOriginalPartFile;

        public void setShouldCallAbortOnCompletion(boolean shouldCallAbortOnCompletion) {
            this.shouldCallAbortOnCompletion = shouldCallAbortOnCompletion;
        }

        public void setIsAnOriginalPartFile(boolean isAnOriginalPartFile) {
            this.isAnOriginalPartFile = isAnOriginalPartFile;
        }

        public int getPartNumber() {
            return this.partNumber;
        }

        public File getPartFile() {
            return this.partFile;
        }

        public long getStart() {
            return this.start;
        }

        public MultipartUploadCallable(int partNumber, File partFile, long start) {
            this.partNumber = partNumber;
            this.partFile = partFile;
            this.start = start;
        }

        @Override
        public UploadPartResult call() throws Exception {
            UploadPartResult result;
            BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(this.partFile));
            inputStream.skip(this.start);
            UploadPartRequest request = (UploadPartRequest)new UploadPartRequest().withBucketName(MultipartUploadManager.this.bucketName).withKey(MultipartUploadManager.this.key).withPartNumber(this.partNumber).withInputStream((InputStream)inputStream).withUploadId(MultipartUploadManager.this.uploadId).withPartSize(this.partFile.length()).withGeneralProgressListener(MultipartUploadManager.this.progressListener);
            boolean doNotDeletePartFile = false;
            try {
                MultipartUploadManager.LOG.info("[PATCHNOTE] uploadPart " + this.partFile.getPath() + " " + this.partFile.length());
                result = MultipartUploadManager.this.s3.uploadPart(request);
                if (MultipartUploadManager.this.closed && this.shouldCallAbortOnCompletion) {
                    MultipartUploadManager.this.s3.abortMultipartUpload(new AbortMultipartUploadRequest(MultipartUploadManager.this.bucketName, MultipartUploadManager.this.key, MultipartUploadManager.this.uploadId));
                }
            }
            catch (AbortedException e) {
                throw e;
            }
            catch (Exception e) {
                MultipartUploadManager.LOG.info("[PATCHNOTE] uploadPart error " + e);
                doNotDeletePartFile = true;
                throw e;
            }
            finally {
                if (!(this.isAnOriginalPartFile && doNotDeletePartFile)) {
                    this.partFile.delete();
                } else {
                    MultipartUploadManager.LOG.info("[PATCHNOTE] Deletion of uploadPart " + this.partFile.getPath() + " averted!");
                }
            }
            return result;
        }
    }

    private class MultipartUploadFutureCallBack
    implements FutureCallback<UploadPartResult> {
        private ListenableFuture<UploadPartResult> future;
        private int partNum;

        public MultipartUploadFutureCallBack(ListenableFuture<UploadPartResult> future, int partNum) {
            this.future = future;
            this.partNum = partNum;
        }

        public void onFailure(Throwable arg0) {
            if (this.future.isCancelled()) {
                MultipartUploadManager.LOG.debug("Multipart Upload for part: " + this.partNum + " cancelled");
            }
        }

        public void onSuccess(UploadPartResult arg0) {
            DateTime endTime = DateTime.now();
            Object object = MultipartUploadManager.this.partNumFutureMapHandle;
            synchronized (object) {
                try {
                    while (MultipartUploadManager.this.partNumFutureMap.get(this.partNum) == null) {
                        MultipartUploadManager.this.partNumFutureMapHandle.wait();
                    }
                }
                catch (InterruptedException e) {
                    throw new RuntimeException("Thread interrupted in multipart upload future callback's onSuccess", e);
                }
            }
            List<MultiPartUploadFuture> multiPartUploadFutures = MultipartUploadManager.this.partNumFutureMap.get(this.partNum);
            MultipartUploadManager.LOG.debug("Total spawned multipart upload futures for partNum: " + this.partNum + " are: " + multiPartUploadFutures.size());
            ArrayList<MultiPartUploadFuture> newMultipartUploadFutures = null;
            for (MultiPartUploadFuture multiPartUploadFuture : multiPartUploadFutures) {
                if (multiPartUploadFuture.getFuture().isDone()) {
                    multiPartUploadFuture.setEndTime(endTime);
                    newMultipartUploadFutures = new ArrayList<MultiPartUploadFuture>();
                    newMultipartUploadFutures.add(multiPartUploadFuture);
                    continue;
                }
                MultipartUploadManager.LOG.debug("Cancelling future for partNum: " + this.partNum + " running for: " + (endTime.getMillis() - multiPartUploadFuture.getStartTime().getMillis()) / 1000 + " s");
                multiPartUploadFuture.getFuture().cancel(true);
                multiPartUploadFuture.getMultiPartUploadCallable().setShouldCallAbortOnCompletion(true);
            }
            if (newMultipartUploadFutures != null) {
                MultipartUploadManager.this.partNumFutureMap.put(this.partNum, newMultipartUploadFutures);
            }
            MultipartUploadManager.LOG.debug("Going to remove " + this.partNum + " from the incomplete part num set");
            MultipartUploadManager.this.incompletePartNums.remove(this.partNum);
        }
    }

    private class MultiPartUploadFuture {
        private long partSize;
        private DateTime startTime;
        private DateTime endTime;
        private Future<UploadPartResult> future;
        private MultipartUploadCallable multipartUploadCallable;

        public long getPartSize() {
            return this.partSize;
        }

        public DateTime getStartTime() {
            return this.startTime;
        }

        public void setStartTime(DateTime startTime) {
            this.startTime = startTime;
        }

        public DateTime getEndTime() {
            return this.endTime;
        }

        public void setEndTime(DateTime endTime) {
            this.endTime = endTime;
        }

        public Future<UploadPartResult> getFuture() {
            return this.future;
        }

        public MultiPartUploadFuture(long partSize, Future<UploadPartResult> future, MultipartUploadCallable multiPartUploadCallable) {
            this.partSize = partSize;
            this.future = future;
            this.multipartUploadCallable = multiPartUploadCallable;
        }

        public MultipartUploadCallable getMultiPartUploadCallable() {
            return this.multipartUploadCallable;
        }
    }

    public static class Builder {
        private MultipartUploadManager multipartUpload = new MultipartUploadManager();

        public Builder withServerSideEncryptionAlgorithm(String s) {
            this.multipartUpload.serverSideEncryptionAlgorithm = s;
            return this;
        }

        public Builder withUploadId(String s) {
            this.multipartUpload.uploadId = s;
            return this;
        }

        public Builder withProgressListener(ProgressListener progressListener) {
            this.multipartUpload.progressListener = progressListener;
            return this;
        }

        public Builder withS3(AmazonS3 s3) {
            this.multipartUpload.s3 = s3;
            return this;
        }

        public Builder withBucketName(String bucketName) {
            this.multipartUpload.bucketName = bucketName;
            return this;
        }

        public Builder withKey(String key) {
            this.multipartUpload.key = key;
            return this;
        }

        public Builder withExecutorService(ListeningExecutorService executorService) {
            this.multipartUpload.executorService = executorService;
            return this;
        }

        public Builder withConf(Configuration conf) {
            this.multipartUpload.conf = conf;
            this.multipartUpload.initializeFromConf(conf);
            return this;
        }

        public Builder withPartSize(long partSize) {
            this.multipartUpload.partSize = partSize;
            return this;
        }

        public MultipartUploadManager build() {
            MultipartUploadManager newMultipartUpload = new MultipartUploadManager();
            newMultipartUpload.serverSideEncryptionAlgorithm = this.multipartUpload.serverSideEncryptionAlgorithm;
            newMultipartUpload.uploadId = this.multipartUpload.uploadId;
            newMultipartUpload.progressListener = this.multipartUpload.progressListener;
            newMultipartUpload.s3 = this.multipartUpload.s3;
            newMultipartUpload.bucketName = this.multipartUpload.bucketName;
            newMultipartUpload.key = this.multipartUpload.key;
            newMultipartUpload.executorService = this.multipartUpload.executorService;
            newMultipartUpload.thFractionPartsCompleted = this.multipartUpload.thFractionPartsCompleted;
            newMultipartUpload.fractionPartAvgCompletionTime = this.multipartUpload.fractionPartAvgCompletionTime;
            newMultipartUpload.partAttempts = this.multipartUpload.partAttempts;
            newMultipartUpload.partSize = this.multipartUpload.partSize;
            newMultipartUpload.conf = this.multipartUpload.conf;
            return newMultipartUpload;
        }
    }

}
EOF

javac -classpath $(find ~/ -name \*.jar | tr '\n' ':') MultipartUploadManager.java
mkdir -p com/amazon/ws/emr/hadoop/fs/s3/
mv *.class com/amazon/ws/emr/hadoop/fs/s3/
jar uf /usr/share/aws/emr/emrfs/lib/emrfs-1.3.0.jar com/amazon/ws/emr/hadoop/fs/s3/*.class
jar uf /home/hadoop/.versions/hbase-0.94.18/lib/emrfs-1.4.0-20150217.191957-3.jar com/amazon/ws/emr/hadoop/fs/s3/*.class
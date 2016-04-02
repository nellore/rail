package edu.jhu.cs;

import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.mapred.TextOutputFormat;
import com.hadoop.compression.lzo.LzopCodec;
import org.apache.hadoop.util.Progressable;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.FileOutputFormat;
import org.apache.hadoop.mapred.RecordWriter;

import java.io.IOException;
import java.io.DataOutputStream;

public class IndexedLzoTextOutputFormat<K, V> extends TextOutputFormat<K, V> {

  @Override
  public RecordWriter<K, V> getRecordWriter(FileSystem ignored,
                                                  JobConf job,
                                                  String name,
                                                  Progressable progress)
    throws IOException {
      String keyValueSeparator = job.get("mapreduce.output.textoutputformat.separator", "\t");
      LzopCodec codec = new LzopCodec();
      codec.setConf(job);

      Path file = FileOutputFormat.getTaskOutputPath(job, name + codec.getDefaultExtension());
      Path indexFile = new Path(file.toString() + ".index");
      FileSystem fs = file.getFileSystem(job);
      FSDataOutputStream fileOut = fs.create(file, progress);
      FSDataOutputStream indexFileOut = fs.create(indexFile, true);
      return new LineRecordWriter<K, V>(new DataOutputStream (codec.createIndexedOutputStream(fileOut, indexFileOut)), keyValueSeparator);
    }

  }

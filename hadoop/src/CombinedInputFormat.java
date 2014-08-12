// Taken verbatim from http://stackoverflow.com/questions/14270317/implementation-for-combinefileinputformat-hadoop-0-20-205
// Also at http://www.bigdataspeak.com/2013/01/how-to-combine-input-files-to-get-to_11.html

package edu.jhu.cs;

import java.io.IOException;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.FileSplit;
import org.apache.hadoop.mapred.InputSplit;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.LineRecordReader;
import org.apache.hadoop.mapred.RecordReader;
import org.apache.hadoop.mapred.Reporter;
import org.apache.hadoop.mapred.lib.CombineFileInputFormat;
import org.apache.hadoop.mapred.lib.CombineFileRecordReader;
import org.apache.hadoop.mapred.lib.CombineFileSplit;

@SuppressWarnings("deprecation")
public class CombinedInputFormat extends CombineFileInputFormat<LongWritable, Text> {

    @SuppressWarnings({ "unchecked", "rawtypes" })
    @Override
    public RecordReader<LongWritable, Text> getRecordReader(InputSplit split, JobConf conf, Reporter reporter) throws IOException {

        return new CombineFileRecordReader(conf, (CombineFileSplit) split, reporter, (Class) myCombineFileRecordReader.class);
    }

    public static class myCombineFileRecordReader implements RecordReader<LongWritable, Text> {
        private final LineRecordReader linerecord;

        public myCombineFileRecordReader(CombineFileSplit split, Configuration conf, Reporter reporter, Integer index) throws IOException {
            FileSplit filesplit = new FileSplit(split.getPath(index), split.getOffset(index), split.getLength(index), split.getLocations());
            linerecord = new LineRecordReader(conf, filesplit);
        }

        @Override
        public void close() throws IOException {
            linerecord.close();
        }

        @Override
        public LongWritable createKey() {
            return linerecord.createKey();
        }

        @Override
        public Text createValue() {
            return linerecord.createValue();
        }

        @Override
        public long getPos() throws IOException {
            return linerecord.getPos();
        }

        @Override
        public float getProgress() throws IOException {
            return linerecord.getProgress();
        }

        @Override
        public boolean next(LongWritable key, Text value) throws IOException {
            return linerecord.next(key, value);
        }

    }
}
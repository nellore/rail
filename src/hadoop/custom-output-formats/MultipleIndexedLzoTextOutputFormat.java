package edu.jhu.cs;

import org.apache.hadoop.mapred.lib.MultipleTextOutputFormat;
import org.apache.hadoop.mapred.RecordWriter;

import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.util.Progressable;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.RecordWriter;

import java.io.IOException;

public class MultipleIndexedLzoTextOutputFormat<K, V> extends MultipleTextOutputFormat<K, V> {

  private IndexedLzoTextOutputFormat<K, V> theIndexedLzoTextOutputFormat = null;

  @Override
  protected String generateFileNameForKeyValue(K key, V value, String name) {
    String k = key.toString();
    return k + "/" + name;
  }

  @Override
  protected K generateActualKey(K key, V value) {
    return null;
  }

  @Override
  protected RecordWriter<K, V> getBaseRecordWriter(FileSystem fs, JobConf job,
      String name, Progressable arg3) throws IOException {
    if (theIndexedLzoTextOutputFormat == null) {
      theIndexedLzoTextOutputFormat = new IndexedLzoTextOutputFormat<K, V>();
    }
    return theIndexedLzoTextOutputFormat.getRecordWriter(fs, job, name, arg3);
  }

}

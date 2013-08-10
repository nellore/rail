package org.myorg;

// import java.io.DataOutputStream;
// import java.io.IOException;
// import java.io.UnsupportedEncodingException;
 
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.mapred.FileOutputFormat;
import org.apache.hadoop.mapred.lib.MultipleTextOutputFormat;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.RecordWriter;
import org.apache.hadoop.mapred.Reporter;
import org.apache.hadoop.util.*;
 
public class MultipleOutputFormat<K, V>
    extends MultipleTextOutputFormat<K, V> {
 
  @Override
  protected String generateFileNameForKeyValue
      (K key, V value, String name) {
      String k = key.toString();
      return k +"/"+name;
  }
}
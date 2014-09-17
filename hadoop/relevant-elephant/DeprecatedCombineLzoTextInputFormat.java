package com.twitter.elephantbird.mapred.input;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;

import com.twitter.elephantbird.mapreduce.input.CombineLzoTextInputFormat;

/**
 * mapred version of {@link CombineLzoTextInputFormat}.
 */
public class DeprecatedCombineLzoTextInputFormat extends DeprecatedFileInputFormatWrapper<LongWritable, Text> {
  public DeprecatedCombineLzoTextInputFormat() {
    super(new CombineLzoTextInputFormat());
  }
}

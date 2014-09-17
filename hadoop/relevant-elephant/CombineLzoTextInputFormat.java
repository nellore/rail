package com.twitter.elephantbird.mapreduce.input;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;

import com.twitter.elephantbird.mapreduce.input.combine.DelegateCombineFileInputFormat;
import com.twitter.elephantbird.mapreduce.input.LzoTextInputFormat;

/**
 * "combine" version of {@link LzoTextInputFormat}.
 */
public class CombineLzoTextInputFormat extends DelegateCombineFileInputFormat<LongWritable, Text> {
  public CombineLzoTextInputFormat() {
    super(new LzoTextInputFormat());
  }
}

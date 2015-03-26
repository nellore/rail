package edu.jhu.cs;

import java.io.UnsupportedEncodingException;
import java.lang.Math;

import org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner;

public class ModPartitioner<K2, V2> extends KeyFieldBasedPartitioner<K2, V2> {

  @Override
  protected int hashCode(byte[] b, int start, int end, int currentHash) {
    try {
      // Assume numbers to try and make the hash a sum of fields
      // This is terrible if there's more than one field, so just make sure there isn't!
      return currentHash + Math.abs(Integer.parseInt(new String(b, start, end - start + 1, "UTF-8")));
    } catch (NumberFormatException e) {
      // Well if that didn't work, return what we would have been returned.
      return super.hashCode(b, start, end, currentHash);
    } catch (UnsupportedEncodingException e) {
      // Forced to handle this
      throw new RuntimeException("The current system does not support UTF-8 encoding!", e);
    }
  }

}

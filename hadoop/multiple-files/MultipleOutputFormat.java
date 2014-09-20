package edu.jhu.cs;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.FileOutputFormat;
import org.apache.hadoop.mapred.lib.MultipleTextOutputFormat;
import org.apache.hadoop.util.*;

public class MultipleOutputFormat<K, V> extends MultipleTextOutputFormat<K, V> {
	
	@Override
	protected String generateFileNameForKeyValue(K key, V value, String name) {
		String k = key.toString();
		return k + "/" + name;
	}
	
	@Override
	protected K generateActualKey(K key, V value) {
		return null;
	}
}

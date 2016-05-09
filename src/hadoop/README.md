Customizations of Hadoop and ElephantBird

* Indexed LZO outputs (without seeking on S3)
* ModPartitioner -- better distribution of per-chromosome work

These are copied into well known directories on the EMR cluster and compiled into Hadoop jars in a bootstrap action.  Hadoop is configured to include them in the appropriate classpath.


Rail developer documentation
============================

This directory is executable, and is the main entry point for Rail-*.

## Philosophy of Rail

Make it easy to run a consistent analysis across many datasets at once.

* User constructs manifest specifying URLs and how they are related biologically
* Rail provides default configuration files for analysis types: e.g., RNA, ChIP, DNA, ...
* Entire analysis should fit in two commands
* Easy to run locally or in various distributed settings
* Output is summarized and organized to facilitate cross-sample queries
* The steps executed are the minimal steps required to output requested deliverables
* Rail is the "entry point" -- Rail-RNA and other pipelines are not directly executable

## Components of Rail

* Dooplicity: parallel abstraction layer
    * Input: JSON, Intermediate/Output: MapReduce outputs
    * Not concerned with specific steps or pipelines
* Rail: genomics layer
    * Input: manifest & config files
    * Understands steps
    * Ensures all tools are installed & executable
    * Uses dooplicity to run pipeline
* Rail pipelines: assay-specific layer
    * Understands pipelines: how to compose steps, how to prune according to deliverables requested

## Current `src` directory structure

* `cloudformation` -- Rail-* CloudFormation templates
* `dooplicity` -- Rail-* parallel abstraction layer
* `hadoop` -- Rail-* customizations of Hadoop, ElephantBird
* `rna` -- Rail-RNA, though with Rail-core mixed in

## Future directory structure

* `src` -- Executable directory; entry point for Rail
* `src/dooplicity` -- Dooplicity!
* `src/dooplicity/tests` -- End-to-end tests for Dooplicity
* `src/rail/core` -- Rail 
    * Generic aspects of driver & installer
* `src/rail/steps` -- scripts to drive individual map or reduce steps in Rail-*
    * Subdirectories are logical groupings of steps (stretches?): `src/steps/coverage`, `src/steps/align`, `src/steps/splicing`
    * A stretch will tend to consist of steps that run consecutively and can be tested as a unit
    * Steps & stretches may be included in many different lines.  E.g. `src/steps/coverage` might be used in both RNA and ChIP lines
* `src/rail/lines`
    * Assay-specific workflows, e.g. `src/rail/lines/rna`
* `src/rail/tests` -- end-to-end tests for Rail

## Testing Dooplicity

Has some unit tests.

Want to add tests that, given a JSON description of a Dooplicity app, checks that results are the same when app is run in all three modes.

Some of these tests will require an IPython-parallel-compatible cluster.

## Testing Rail

### Unit tests

In theory, this is how you run all the unit tests:

```
for i in `find . -name '*.py'` ; do python $i --test ; done
```

### Stretch tests

Could be disguised as unit tests so that above command also runs stretch tests.  E.g. `src/rail/steps/coverage/test_coverage.py` with tests for the `coverage` stretch and no other functionality.

### End-to-end tests

The tutorial is the closest thing to an end-to-end test.

* See `rail/ex` and docs.rail.bio repo

But there's no real check being done on the output.

## Counters

Counters are an important debugging asset.  Currently, Dooplicity does not recognize or compile counters.  Suggestions for adding counters:

* Counters and logging and keepalive all tied together?
* Counters periodically flushed in coordination with keepalive
* Counters finally flushed using Python atexit
* Up-to-date thread-local counters maintained in an object that "gets passed around" so that unit tests can interrogate how a specific function updated the counters 

## TODO

* RNA tri-mode line test
* Dooplicity tri-mode line test
* Refactor to enable Rail-*
* Dooplicity unit tests
* Rail step unit tests
* Rail stretch tests
* Logger/Counters
* Possibly relegate `eval` to separate repo

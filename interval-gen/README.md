# Generating Intervals

This folder illustrates the generation of reduced rounding intervals
for various elementary functions. To generate the intervals for any
new elementary function, one needs to add the definition of the following functions.

```
GuessInitialLbUb
Range Reduction
OutputCompensation
compute_special_case
```

## Building the intervals for Log2

To generate the intervals for the elementary function Log2, build the
interval generator as follows.

```
make Log2
```

## Executing the Interval Generator for Log2

To execute the interval generator, provide the name of the interval
file and the path to the oracle file with the round-to-odd result for
all inputs.

```
./Log2 log2-intervals.txt <log2_oracle_file>
```



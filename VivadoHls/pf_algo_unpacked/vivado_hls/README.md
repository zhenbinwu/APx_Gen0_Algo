Need to place GlobalCorrelator_HLS code (or a link) in src. 

Eg.
```
cd src
git clone https://github.com/drankincms/GlobalCorrelator_HLS.git -b dev
```

This will use PF+PUPPI as a combined algo block. mp7wrapped_pfalgo3_full() is the top function.

To run:
```
vivado_hls -f run_hls.tcl
```

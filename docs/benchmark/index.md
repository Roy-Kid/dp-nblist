# Performance and Benchmark

In this section, we will present the performance of the library and compare it with other libraries. The performence covers speed and memory usage. The results will be organized in the following way:

* Scaling:
    It will test the performance of the library with different size of data. The size of data will be increased by exponential of 2 and the data will be generated randomly with consistent density 0.6

* Structure:
    It will test the performance for different structrues.
    * Interfacial: solid-liquid system;
    * Protein: protein in water;
    * Heterogeneous: system with different shape;

* Latency:
    It will test speed for re-build from scratch. It is important for small molecule training and docking.
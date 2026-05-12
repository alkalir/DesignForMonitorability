# MONICATool

This folder contains the resources related to the classification of on-chip monitoring systems (OCMSs) using the [MONICA tool](https://monicatool.cloud/).

## Folder Structure

- **parametersDescription:** contains the description of all the parameters used to classify OCMSs. For each parameter, the description includes the definition, the set of possible values, and, where applicable, examples.
- **underClassification:** contains the OCMSs that are currently under classification.

## How to Classify an OCMS

1. Pick one OCMS from the `underClassification` folder.
2. Read the description of each parameter in the `parametersDescription` folder.
3. Create a text file named after the OCMS, reporting the value assigned to each parameter. Use the following format — if a parameter has multiple values, list them separated by a comma:
```
Purpose: Performance, Security
Hardware Target: RISC-V
Policy: Non-Blocking
...
```
4. Submit the classification by sending the file to [giacomo.valente@univaq.it](mailto:giacomo.valente@univaq.it).

## Suggested Order for Filling in the Parameters

We suggest filling in the parameters in the following order, going from a high-level characterization of the monitoring system down to its costs and implementation details:

1. Purpose — what the monitor is designed for
2. Hardware Target — the hardware platform being monitored
3. Software Platform — the software layer where the monitored application runs
4. Organization — the MB-MI-MP structure of the monitoring system
5. Technique — how the monitoring system collects data
6. Metrics — what the monitoring system measures
7. Granularity — the level of detail of the collected events
8. Policy — whether the monitor is blocking or non-blocking
9. Multi-Thread — support for multi-threaded workloads
10. Multi-Core — support for multi-core processors
11. System-Wide — support for heterogeneous GPP+SPP environments
12. Synchronization — how synchronization among MB, MI, and MP is achieved
13. Resource Sharing — shared resources among MB, MI, and MP
14. Exceptions — whether the monitor can raise exceptions/interrupts
15. Record-Replay — whether the monitor supports record and replay
16. Detection Latency — the time between event occurrence and monitoring information availability
17. Performance Degradation — the timing overhead introduced by the monitor
18. Memory Space — the memory occupation of the monitor
19. Physical Cost — the area occupied by the monitor
20. Power Dissipation Cost — the power consumption increase due to the monitor
21. Design Cost — the effort needed to design and implement the monitor
22. Specification Modification — how much the monitored application needs to be modified
23. Accuracy — the correctness of the measurements or classifications
24. Portability — the possibility to reuse the monitor on different targets
25. Extensibility — the capability to extend the monitor with new functionalities
26. Availability — whether the monitor implementation is available
27. Type of Physical Implementation — the type of physical platform where the monitor runs
28. Automatic Generation — whether the monitor can be automatically generated
29. Runtime Management — whether the monitor offers runtime management capabilities

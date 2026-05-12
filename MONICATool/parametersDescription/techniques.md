# Technique

The **Technique** parameter refers to the technique employed to perform a monitoring action.

This parameter can be set with one or more values from the list below.

## Possible Values

- **Source Level Code Instrumentation:** this type of instrumentation is related to the addition of source code directives within the source code of the monitored application, in order to extract event instances.
- **Binary Level Code Instrumentation:** this type of instrumentation is related to the addition of directives at binary level, in order to extract event instances. This can be statically performed (Static Binary Level Code Instrumentation) or dynamically performed (Dynamic Binary Level Code Instrumentation).
- **Sampling:** a monitoring system that adopts the sampling approach wraps itself around the execution of an application, taking control of the program flow, and then pausing the execution at specific points to record the current state of the system (such as reading the program counter register).
- **Hardware Performance Counters:** they represent hardware elements composed of a counting register with a control system, able to collect some metrics related to performance events.
- **Performance Monitoring Units:** they are dedicated hardware elements built inside a processor to measure its performance parameters. The difference with a hardware performance counter is that a performance monitoring unit is not necessarily based on events counting operations.
- **Hardware Trace Buffers:** they are dedicated hardware elements able to trace an application execution, and either store those traces or output them.
- **Custom:** the system uses a custom technique to perform the monitoring action.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

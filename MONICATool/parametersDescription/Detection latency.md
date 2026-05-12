# Detection Latency

The **Detection Latency** parameter refers to the elapsed time between the instant in which the event instance to be observed happens and the instant where the monitoring information is available.

This parameter can be set with an absolute value. 

## Possible Values

- **Maximum Detection Latency in clock cycles (cc):** a number or a range related to the detection latency value expressed in clock cycles. If the detection latency is provided in ms, it must be scaled to clock cycles by using the frequency of operation.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- 4 cc (e.g., a hardware monitoring system)
- 100 cc (e.g., a software monitoring system running on bare-metal)

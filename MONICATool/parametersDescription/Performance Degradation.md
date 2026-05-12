# Performance Degradation

The **Performance Degradation** parameter refers to the maximum timing overhead caused by the use of the monitoring system.

This parameter can be set with an absolute value.

## Possible Values

- **Maximum execution time increase (cc):** the maximum increase in execution time, expressed in clock cycles, due to the use of the on-chip monitoring system.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- 10 cc (e.g., a hardware monitoring system with minimal interference)
- 500 cc (e.g., a software monitoring system on a bare-metal application)

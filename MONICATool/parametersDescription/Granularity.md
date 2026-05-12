# Granularity

The **Granularity** parameter refers to the nature of the events that can be collected by the monitoring system.

This parameter can be set with one or more values from the list below.

## Possible Values

- **Architectural:** events observable from the software level, e.g., number of retired instructions, branches, cycles.
- **Micro-Architectural:** events observable only at the hardware level, e.g., cache accesses, branch prediction, TLB accesses.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

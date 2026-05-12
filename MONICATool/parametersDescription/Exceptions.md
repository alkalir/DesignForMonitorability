# Exceptions

The **Exceptions** parameter refers to the capability of the on-chip monitoring system to raise an exception/interrupt.

This parameter can be set with one value from the list below.

## Possible Values

- **YES:** the monitoring system is capable of raising an exception/interrupt.
- **NO:** the monitoring system is not capable of raising an exception/interrupt.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Background

In digital computers, an interrupt is a signal sent to a processor to require its attention. An interrupt alerts the processor and serves as a request for the processor to interrupt the currently executing code when permitted, so that the event can be processed in a timely manner. If the request is accepted, the processor generally responds by suspending its current activities, saving its state, and executing a function called interrupt handler (or an interrupt service routine, ISR) to deal with the event.

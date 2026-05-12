# Type of Physical Implementation

The **Type of Physical Implementation** (TPI) parameter refers to the physical implementation of the monitoring system. In other words, TPI refers to the type of physical platform where the monitoring system can perform the monitoring action.

This parameter can be set with one or more values from the list below.

## Possible Values

- **Completely Fixed (CF):** this TPI refers to hardwired implementations. They can be, for example, a hardwired General Purpose Processor (GPP, such as Intel ATOM) or a hardwired Single Purposed Processor (SPP, such as a GPU) without any available reconfigurable logic.
- **Reconfigurable (RCF):** this TPI refers to reconfigurable logic based implementations, such as implementations on Field Programmable Gate Arrays (FPGAs).
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- CF (e.g., a monitoring system running on an Intel ATOM processor)
- RCF (e.g., a monitoring system implemented on an FPGA)

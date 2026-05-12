# Software Platform

The **Software Platform** parameter refers to on-chip monitoring systems that target a software monitored behaviour (MB = SW, with respect to the monitoring reference model, see the [Organization](organization.md) parameter for more details), and it refers to the software layer where the application to be monitored is executed.

This parameter can be set with one value from the list below.

## Possible Values

- **Bare-metal:** the monitored application runs directly on the hardware without any operating system.
- **Operating-system:** the monitored application runs on top of an operating system.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

# Policy

The **Policy** parameter refers to whether the on-chip monitoring system, while applying its monitoring process, is blocking or non-blocking for the execution of the monitored application.

This parameter can be set with one value from the list below.

## Possible Values

- **Blocking:** the monitoring system blocks the execution of the monitored application during the monitoring process.
- **Non-Blocking:** the monitoring system does not block the execution of the monitored application during the monitoring process.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

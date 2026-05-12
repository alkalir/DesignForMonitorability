# Specification Modification

The **Specification Modification** parameter refers to how much the behaviour of the monitored application needs to be modified in order to adopt the monitoring system. The value represents the percentage of the total system specification that needs to be modified, ranging from 0% (no modifications needed) to 100% (the specification shall be completely rewritten).

This parameter can be set with an absolute value.

## Possible Values

- **value (%):** the percentage of the total system specification that needs to be modified in order to adopt the monitoring system, from 0% to 100%.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- 0% (e.g., a hardware monitoring system that does not require any modification to the monitored application)
- 100% (e.g., a monitoring system that requires a complete rewrite of the monitored application specification)

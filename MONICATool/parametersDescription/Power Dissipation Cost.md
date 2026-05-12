# Power Dissipation Cost

The **Power Dissipation Cost** parameter refers to the maximum power dissipation cost related to the use of the monitoring system.

This parameter can be set with an absolute value.

## Possible Values

- **Maximum power dissipation increase (mW):** the maximum increase in power dissipation, expressed in milliwatts, due to the use of the on-chip monitoring system.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- 5 mW (e.g., a lightweight hardware monitor)
- 50 mW (e.g., a software monitoring system)

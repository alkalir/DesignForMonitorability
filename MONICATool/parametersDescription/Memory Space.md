# Memory Space

The **Memory Space** parameter refers to the data memory requirements for the storage of raw information and the storage of monitoring information. Moreover, it also includes the instruction and data memory requirements for the storage of necessary information for the correct operations of the monitoring system.

This parameter can be set with one or both of the following values.

## Possible Values

- **Data Memory occupation (Bytes):** the maximum value of data memory occupied by the monitoring system, expressed in Bytes.
- **Instruction Memory occupation (Bytes):** the maximum value of instruction memory occupied by the monitoring system, expressed in Bytes.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- Data Memory occupation: 1024 Bytes (e.g., a lightweight hardware monitor)
- Instruction Memory occupation: 4096 Bytes (e.g., a software monitoring system)

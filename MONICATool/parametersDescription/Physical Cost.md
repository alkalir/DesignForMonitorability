# Physical Cost

The **Physical Cost** parameter refers to the maximum cost, in terms of occupied area, of the physical components added to the platform due to the implementation of the on-chip monitoring system.

This parameter can be set with one or more values from the list below.

## Possible Values

- **Maximum occupied area in terms of Flip-Flops:** the maximum number of Flip-Flops occupied by the on-chip monitoring system.
- **Maximum occupied area in terms of LUTs:** the maximum number of Look-Up Tables (LUTs) occupied by the on-chip monitoring system.
- **Maximum occupied area in terms of BRAMs:** the maximum number of Block RAMs (BRAMs) occupied by the on-chip monitoring system.
- **Maximum occupied area in terms of DSPs:** the maximum number of Digital Signal Processors (DSPs) occupied by the on-chip monitoring system.
- **Maximum occupied area in terms of mm²:** the maximum silicon area in mm² occupied by the on-chip monitoring system.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Examples

- Maximum occupied area in terms of LUTs: 128 (e.g., a lightweight hardware monitor on FPGA)
- Maximum occupied area in terms of mm²: 0.05 mm² (e.g., a monitor implemented on an ASIC)

# Hardware Target

The **Hardware Target** parameter refers to the hardware target monitored by the on-chip monitoring system.

This parameter can be set with one or more values from the list below. General-purpose Processor and Single-purpose Processor are generic categories to be used when the precise architecture is not known or not specified. ARM Cortex A9, Microblaze, and RISC-V are specific GPP architectures. AMBA3 AHB-Lite refers to a bus architecture.

## Possible Values

- **General-purpose Processor (GPP):** a programmable processor designed to execute programs expressed with an instruction set. To be used when the precise architecture is not specified.
- **Single-purpose Processor (SPP):** a processor designed to execute a specific task, such as an accelerator. To be used when the precise architecture is not specified.
- **ARM Cortex A9:** a specific GPP with ARM Cortex A9 architecture.
- **Microblaze:** a specific GPP with Xilinx MicroBlaze architecture.
- **RISC-V:** a specific GPP with RISC-V architecture.
- **AMBA3 AHB-Lite:** a bus architecture based on the AMBA3 AHB-Lite standard.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

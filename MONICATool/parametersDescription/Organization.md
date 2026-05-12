# Organization

The **Organization** parameter refers to the relation, in terms of implementation strategy, between the monitored application and the on-chip monitoring system. The relation is expressed through the monitoring reference model.

The monitoring reference model abstracts the generic monitoring process with three parts: a Monitored Behaviour (MB), representing the application under monitoring, and a Monitoring Infrastructure (MI) and a Monitoring Processor (MP), both part of the on-chip monitoring system. They are detailed in the following:

- **Monitored Behaviour (MB):** represents the tasks, executing on a target system, that have to be monitored with an on-chip monitoring system. MB can be implemented in hardware (executed by dedicated and not reprogrammable hardware architecture) or software (executed by a programmable processor).

- **Monitoring Infrastructure (MI):** represents all the necessary mechanisms, part of an on-chip monitoring system, able to extract raw information. A software mechanism uses the resources of the target under analysis to perform the raw information extraction. A hardware mechanism has a dedicated hardware architecture to perform the raw information extraction.

- **Monitoring Processor (MP):** represents all the necessary tasks, part of an on-chip monitoring system, that by using raw information as inputs, apply some algorithms to provide monitoring information organized in metrics. MP can be implemented in hardware or software.

<div align="center">

![organization](https://i.imgur.com/nJbR3yb.png)

</div>

This parameter can be set with one or more values from the list below, expressed as a triple (MB, MI, MP) indicating the type of implementation of each part.

## Possible Values

- **SW-SW-SW**
- **SW-SW-HW**
- **SW-HW-SW**
- **SW-HW-HW**
- **HW-SW-SW**
- **HW-SW-HW**
- **HW-HW-SW**
- **HW-HW-HW**
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

## Example

- SW-HW-SW: the MB is a software running on a programmable processor, MI is implemented with dedicated hardware elements, and MP is again a software. This triple can describe a hardware monitoring system that extracts information about a software workload execution without introducing software overhead, and that aggregates this information, for example to get an estimation of worst-case execution time, using a software executing on the same target.

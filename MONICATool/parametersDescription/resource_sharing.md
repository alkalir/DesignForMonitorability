# Resource Sharing

The **Resource Sharing** parameter refers to the elements that are shared between MB, MI and MP (with respect to the monitoring reference model, see the [Organization](organization.md) parameter for more details).

The analysis considers possible sharing among physical and logical resources of MB, MI and MP, as shown in the figure below:

<div align="center">

![resource-sharing](https://i.imgur.com/0JKfVvN.png)

</div>

This parameter can be set with one or more values from the list below.

## Possible Values

- **MB-MI: Processing Resources**
- **MB-MI: Instruction Memory**
- **MB-MI: Data Memory**
- **MB-MI: No Sharing**
- **MB-MP: Processing Resources**
- **MB-MP: Instruction Memory**
- **MB-MP: Data Memory**
- **MB-MP: No Sharing**
- **MI-MP: Processing Resources**
- **MI-MP: Instruction Memory**
- **MI-MP: Data Memory**
- **MI-MP: No Sharing**
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

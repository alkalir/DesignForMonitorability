# Synchronization

The **Synchronization** parameter refers to the required synchronization among elements that implement MI, MB, and MP, in order to correctly perform the monitoring action.

This parameter can be set with one value from the list below.

## Possible Values

- **Automatic:** the synchronization does not require any additional actions, as it is guaranteed by the use of the monitoring system itself. This is the case, for example, of a user space operating system application monitored with the technique of source level code instrumentation, where the compilation process ensures the correct synchronization among MB, MI and MP.
- **Manual:** a hardware/software custom component is required to guarantee the synchronization.
- **Not Applicable:** This value should never be set.
- **Not Declared:** this parameter is not declared for the considered monitoring requirements, i.e., no value has been specified for it. This value must be set in case the other ones are not set.

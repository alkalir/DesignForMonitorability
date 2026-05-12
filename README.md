# Design for Monitorability

> If you use the proposed approach or the MONICA tool in your work, please cite \[1\].

This repository contains the code and resources associated with the **Design for Monitorability** project, whose goal is to develop an approach and associated framework to integrate on-chip monitoring systems into the HW/SW co-design phase.

The HW/SW co-design tool of reference is [Hepsycode](https://hepsycode.github.io/), while the taxonomy and parameter set developed to describe on-chip monitoring systems is available at [MONICATool](https://monicatool.cloud/).

---

## Repository Structure

### FFG — FIR-FIR-GCD Synthetic Benchmark

The `FFG` folder contains the experimental results related to a synthetic benchmark (FIR-FIR-GCD) used to validate the proposed approach, as described in \[1\].

### Pacemaker

The `pacemaker` folder contains the experimental results related to the Pacemaker application, also used for the validation of the proposed approach in \[1\].

### Basic Hardware Elements

The `basic_hardware_elements` folder provides the HDL descriptions of the basic hardware elements considered in the experiments (processors, physical links, and memories), enabling full reproducibility of the tests presented in \[1\].

### MONICATool — Parameters Description

The `MONICATool/parametersDescription` folder contains the descriptions of the parameters used in the [MONICA tool](https://monicatool.cloud/) to characterize on-chip monitoring systems. Each parameter is documented in a dedicated Markdown file.

---

## References

\[1\] Giacomo Valente, Vittoriano Muttillo, Luigi Pomante, Daniele Frigioni, and Tania Di Mascio. 2025. *A New HW/SW Co-Design Approach for Monitored Systems-on-Chip Development*. ACM Trans. Embed. Comput. Syst. 24, 6, Article 172 (November 2025), 37 pages. https://doi.org/10.1145/3769075

\[2\] Tania Di Mascio, Federica Caruso, Laura Tarantino, and Giacomo Valente. 2021. *MONICA Vision: An Approach, a Model and the Interactive Tools for Cyber-Physical Systems Designers*. In Proceedings of the 14th Biannual Conference of the Italian SIGCHI Chapter (CHItaly '21). Association for Computing Machinery, New York, NY, USA, Article 33, 1–5. https://doi.org/10.1145/3464385.3464778

\[3\] Caruso, F., Di Mascio, T., Peretti, S., Pomante, L., Valente, G. 2022. *MONICA "On-the-Job" Technology-Enhanced Learning Environment: An Empirical Evaluation*. In: De la Prieta, F., et al. Methodologies and Intelligent Systems for Technology Enhanced Learning, 11th International Conference. MIS4TEL 2021. Lecture Notes in Networks and Systems, vol 326. Springer, Cham. https://doi.org/10.1007/978-3-030-86618-1_6

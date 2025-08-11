# Portfolio Optimization

Project name: 
* Annealed Portfolio

Team Members:
* Hao Mack Yang (gst-Z2od5sCatjVtRvH)

Presentation:
* [https://studio.youtube.com/video/ZSvtRcaG_Vc/livestreaming](https://youtube.com/live/ZSvtRcaG_Vc?feature=share)

This **Washington Institute for STEM Entrepreneurship and Research** project demonstrates an alternative quantum computing method to optimize the portfolios based on the investor's preferences.

The goal is to minimize an objective function, the sum of the squared difference between the actual amount against the target amount of a target characteristic by risk bucket membership and characteristic while subject to guardrail constraints regarding the actual amounts. The minimization of this objective function allows us to maximize the expected utility of our portfolio. Our output is a set of binary decision variables indicating whether the particular bond is included in the portfolio.

Our approach uses the D-WAVE quantum annealing method, as opposed to the variational quantum algorithm implemented in Qiskit or Xanadu. Quantum annealing is especially well-suited for scalability and handling quadratic optimization problems, even with constraints. Unfortunately, we are not able to get access to the actual quantum processor through LEAP, so we ended up using one of D-WAVE's simulators. Additionally, constrained binary quadratic optimization problem still requires a hybrid approach.

The structure of this submission is such that with elementary knowledge of Python or Jupyter plus the introduction of virtual environments, one can interactively see the process of data gathering, guardrail inputting, problem conversion, and simulated annealing, and presentation with one single click of the play button.

With the exception of the README and presentations, all content that belongs to the team is in the `implementation` folder relative to the root working directory. Inside contains another README.md regarding the mathematical formulation, an `.ipynb` file depicting a large dataset run with randomly set bounds, a slightly modified `.ipynb` file from the `lp` bounds, a `.py` file depicting repeated runs of the small dataset runs, and the `bond_data1.json` generated for listing the selected attributes of the reduced dataset.

Upon running the script, a `dwave_result.csv` file will be generated at the `implementation` working directory (caveat: it is relative to the working directory of VS Code if you are running on VS Code), and various results files depicting the energy value, bond selection, and constraint satisfaction of each entry ordered by the D-WAVE energy rankings (the lowest energy selection is considered the best solution.)

We have supplied the `result` and `result2` images, depicting the energy and secondary constraint satisfaction trends ranked by D-WAVE energy rating, and the `qaoamatrix.png` images, depicting the symmetrical matrix representing the binary QAOA matrix.

Please bear in mind that there are some limitations present in this project, especially concerning the potential differences between the actual boundaries used by the reference method against our method. The structure of the code is designed such that there should be minimal debugging and code redundancy when correcting the code to match the reference method to evaluate against.

Word count: 460
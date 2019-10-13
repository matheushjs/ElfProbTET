# ElfProbTET
Probabilistic Treatment of Execution Times

This is the main repository of my computer science course conclusion project 1.

<div class="container text-justify">
  <h1>Objectives</h1>
  <p>
    The objective of this project focuses on the scheduling problem, but in a way that may also contribute to the two other areas mentioned in the previous section.
    More specifically, focus on the problem of <i>DAG scheduling</i>[1],
      in which tasks are organized as a dependency graph and the objective is to find the best order for their execution on the available computational resources,
      aiming to achieve the lowest overall execution time (also called <i>makespan</i>).
    The problem is made more complicated because, as often happens, the DAG being scheduled will be executed alongside multiple other DAGs (varying with time),
      and these will make use of the same available resources.
  </p>
  <p>
    The scheduling problem has often been approached in a deterministic way, using the average execution time of tasks to perform the scheduling.
    Some probabilistic models have also been proposed for the problem (in this case it is called <i>stochastic scheduling</i>),
      mostly under the assumption that the probability distribution of tasks is known,
      such as done by Li and Antonio[2] and Zheng and Sakellariou[1].
    In both of these papers, in order to validate their results, the authors performed simulations where the the distributions of tasks are somewhat arbitrarily chosen to be either normal or uniform.
    In this project, we would like to further investigate whether these are reasonable choices of distributions for execution times.
    Under this light, we propose the following:
  </p>
  <p>
    <b>Hypothesis 1.</b>
    Given a computer architecture, it is possible to determine a minimal probability model that can be fit to the execution time of any program for a fixed input.
  </p>
  <p>
    In order to investigate this hypothesis, a set of programs will be implemented and will undergo experiments to generate samples of execution times.
    Based on these samples, suitable probability models will be determined, and inference (conventional and/or bayesian) will be performed to find the parameters of these distributions for each sample of execution times.
  </p>
  <p>
    If Hypothesis 1 is true, simulations will be performed to validate or reproduce the results given in existing studies regarding stochastic scheduling.
    If it is not true, the reasons thereof will be investigated, and the hypothesis will be narrowed down to more specific classes of programs and machine circumstances,
      which will undergo the same procedure exposed above.
  </p>
  <p style="margin-top: 2rem;">
    [1] - ZHENG, W.; SAKELLARIOU, R. Stochastic dag scheduling using a monte carlo approach. Journal of Parallel and Distributed Computing, Elsevier, v. 73, n. 12, p. 1673–1689, 2013.
    <br>
    [2] - LI, Y. A.; ANTONIO, J. K. Estimating the execution time distribution for a task graph in a heterogeneous computing system.
      In: IEEE. Proceedings Sixth Heterogeneous Computing Workshop (HCW’97). [S.l.], 1997. p. 172–184.
  </p>
</div>

<div class="container">
  <h1>Results</h1>
  <ul>
    <li><a target="_blank" href="https://cran.r-project.org/web/packages/elfDistr/index.html">R Package "elfDistr"</a></li>
    <li><a target="_blank" href="https://github.com/matheushjs/ElfProbTET">Project GitHub Repository</a></li>
  </ul>
</div>

<div class="container">
  <h1>Funding & Support</h1>
  <ul>
    <li>
      Computational resources from <a target="_blank" href="http://www.cemeai.icmc.usp.br/">CeMEAI</a> (FAPESP grant 2013/07375-0).
    </li>
  </ul>
</div>

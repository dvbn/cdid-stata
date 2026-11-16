# Attribution and licenses

`cdid` for Stata is distributed under version 3 of the GNU General Public License. The complete license is in [LICENSE](LICENSE).

Parts of the code are adapted from [cdid for R](https://github.com/dvbn/cdid), version 1.1.0.0 (commit `1039120668cceabd8111cad82d32ab3927c2b94f`), licensed GPL-3, whose authors are David Benatia, Christophe Bellégo, Joel Cuerrier and Vincent Dortet-Bernadet. The headers of the source files identify the adapted code. The Stata results are compared with version 1.2.0.1 of the R package.

Ordinary mode uses [csdid](https://github.com/pedrohcgs/csdid-stata), version 2.0.2 (commit `be5d272debb1c98e18c7fba79f51b2b3effd3fc4`), by Brantly Callaway, Fernando Rios-Avila and Pedro H. C. Sant'Anna, which is distributed unchanged in the `csdid` folder and installed separately. The graph layout of `cdid_plot` is adapted from `csdid_plot` 2.0.2. Their MIT license, with its copyright notice, is reproduced in [references/csdid-MIT-LICENSE.txt](references/csdid-MIT-LICENSE.txt); it applies to that code and does not replace the license of `cdid`.

No Stata executable, license file or manual is distributed with `cdid`.

## References

Bellégo, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained Difference-in-Differences. *Journal of Econometrics*, 248, 105783. <https://doi.org/10.1016/j.jeconom.2024.105783>

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with Multiple Time Periods. *Journal of Econometrics*, 225(2), 200-230. <https://doi.org/10.1016/j.jeconom.2020.12.001>

Some source comments in this distribution have been edited for publication. Executable statements, numerical conventions, license terms and copyright notices are unchanged.

# Non-convex_Approach_to_Joint_Sensor_Calibration_and_Spectrum_Estimation
Direction of Arrival (DoA) estimation with auto-sensor calibration



A. INTRODUCTION
Blind sensor calibration for spectrum estimation is the problem of estimating the unknown sensor calibration parameters as well as the parameters-of-interest of the impinging signals simultaneously from snapshots of measurements obtained from an array of sensors. In this paper, we consider blind phase and gain calibration (BPGC) problem for direction-of-arrival estimation with multiple snapshots of measurements obtained from an uniform array of sensors, where each sensor is perturbed by an unknown gain and phase parameter. Due to the unknown sensor and signal parameters, BPGC problem is a highly nonlinear problem. Assuming that the sources are un- correlated, the covariance matrix of the measurements in a perfectly calibrated array is a Toeplitz matrix. Leveraging this fact, we first change the nonlinear problem to a linear problem considering certain rank-one positive semidefinite matrix, and then suggest a non-convex optimization approach to find the factor of the rank-one matrix under a unit norm constraint to avoid trivial solutions. Numerical experiments demonstrate that our proposed non-convex optimization approach provides better or competitive recovery performance than existing methods in the literature, without requiring any tuning parameters. 



B. GOALS
Direction of Arrival (DoA) estimation with auto-sensor calibration



C. CONTENTS: matlab codes


D. Prerequisite software (or solver) to run the MATLAB code: Manopt
https://www.manopt.org


E. REFERENCES:
[1] M. Cho, W. Liao and Y. Chi, "A Non-Convex Approach To Joint Sensor Calibration And Spectrum Estimation," in Proceedings of IEEE Statistical Signal Processing Workshop (SSP), 2018, pp. 398-402.
https://ieeexplore.ieee.org/abstract/document/8450691



E. CONTECT INFO.: Myung (Michael) Cho 
- Email: michael.myung.cho@gmail.com
- Homepage: https://sites.google.com/view/myungcho/home

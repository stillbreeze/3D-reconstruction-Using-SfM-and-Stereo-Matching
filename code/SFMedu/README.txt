Run SFMedu2.m

It calls full_pipeline() with several parameters 
including dataset parameter, visualization parameter, visualize
reprojection images boolean, dense reconstruction boolean (off by default to save time)
and runs on all datasets

Can adjust the do_sequential_match parameter to switch between MST matching and sequential matching.

bundleAdjustmentFull adjusts all intrinsics and extrinsics, along with interative motion only & structure only steps
reprojectionErrorFull computes reprojection error with full intrinsic matrix



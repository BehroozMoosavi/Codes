function [yhat_level, yhat_z_raw] = back_transform(yhat_z_std, y_lag1, z_mean, z_std)
%% back_transform.m  |  FILE 12 of 30
%  Un-standardize z_hat and invert the log-change transformation.
%  Step 1: z_raw  = z_mean + z_std * z_std_hat
%  Step 2: y_hat  = y_{t-1} * exp(z_raw / 100)

    yhat_z_raw = z_mean + z_std .* yhat_z_std;
    yhat_level = y_lag1 .* exp(yhat_z_raw / 100);
end

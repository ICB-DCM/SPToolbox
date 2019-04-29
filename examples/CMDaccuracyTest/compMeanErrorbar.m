function [dApp_MC, dApp_MC_neg, dApp_MC_pos] = compMeanErrorbar(model, ...
         xi,estruct, op_SP, PointN, dim_out, dim_in, MCP)
[dApp_MC.my, dApp_MC.Cy, dApp_MC_neg.my, dApp_MC_neg.Cy, dApp_MC_pos.my,...
    dApp_MC_pos.Cy] = deal(zeros(size(PointN, 2), 1));
for i_samples = 1:size(PointN, 2)
    op_SP.n_samples = PointN(i_samples);
    rand_num = 100;
    dApprand_MC.my = zeros(rand_num, dim_out);
    dApprand_MC.Cy = zeros(rand_num, dim_out, dim_out);
%% compute rotate data
    for jrand = 1:rand_num
        switch op_SP.approx
            case 'halton'
                op_SP.skip = (jrand - 1) * op_SP.n_samples + 100;
            otherwise
                op_SP.angle = createRandomRotationMatrix(dim_in);
        end
        Apprand = getSigmaPointApp(model,xi,estruct,op_SP);
        dApprand_MCt = structminus(Apprand, MCP);
        [dApprand_MC.my(jrand, :), dApprand_MC.Cy(jrand, :, :)] = average_time(dApprand_MCt);
    end
    [dApprand_MC] = average_difference(dApprand_MC, dim_out);        
    [dApp_MC.my(i_samples), dApp_MC.Cy(i_samples)] = average_time(dApprand_MC);
    error_App_MC.my = sort(dApprand_MC.my - dApp_MC.my(i_samples));
    error_App_MC.Cy = sort(dApprand_MC.Cy - dApp_MC.Cy(i_samples));
    error_App_MC = structfun(@(x) x(6:95), error_App_MC, 'UniformOutput', false);
    dApp_MC_neg.my(i_samples) = abs(min(error_App_MC.my));
    dApp_MC_neg.Cy(i_samples) = abs(min(error_App_MC.Cy));
    dApp_MC_pos.my(i_samples) = abs(max(error_App_MC.my));
    dApp_MC_pos.Cy(i_samples) = abs(max(error_App_MC.Cy));
end
end
function [status] = Test_Bench (Num_Scenarios,Num_SubScenarios,Num_tar,Num_int,Override)
Num_TotalScenarios = Input_DataSet_Generation(Num_Scenarios,Num_SubScenarios,Num_tar,Num_int,Override);
status = 0;
for i = 1:Num_TotalScenarios
    TestBed_Function(i);
end
status = 1;
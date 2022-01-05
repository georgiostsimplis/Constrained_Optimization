function [x,lambda]=EqualityQPSolver(H,g,A,b,solver)
expectedsolvers = {'NullSpace','RangeSpace','LUdense',...
    'LUsparse','LDLdense','LDLsparse'};

solvername = validatestring(solver,expectedsolvers,6);
switch solvername
    case {'LUdense'}
        [x,lambda]=EqualityQPSolverLUdense(H,g,A,b);
    case {'LUsparse'}
        [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b);
    case {'LDLdense'}
        [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b);
    case {'LDLsparse'}
        [x,lambda]=EqualityQPSolverLDLsparse(H,g,A,b);
    case {'NullSpace'}
        [x,lambda]=EqualityQPSolverNS(H,g,A,b);
    case {'RangeSpace'}
        [x,lambda]=EqualityQPSolverRS(H,g,A,b);
    otherwise
        error('Unknown Solver')
end
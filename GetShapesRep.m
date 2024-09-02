%Define Mai Function
function main()
    ShapeToGenerate = 20;
    ranges = defineRanges();
    generateShapes(ShapeToGenerate, ranges);
end

%Function For The Ranegs 
function ranges = defineRanges()
    ranges.A = [0.1, 10];
    ranges.B = [0.1, 10];
    ranges.M = [4, 4];
    ranges.N = [0.1, 100];
end

%Generate Shpes
function generateShapes(ShapeToGenerate, ranges)
    for i = 1:ShapeToGenerate
        params = generateRandomParameters(ranges);
        createAndDisplayShape(params, i, ShapeToGenerate);
    end
end

%Parameters
function params = generateRandomParameters(ranges)
    params.a1 = randRange(ranges.A);
    params.b1 = randRange(ranges.B);
    params.m1 = randRange(ranges.M);
    params.n11 = randRange(ranges.N);
    params.n21 = randRange(ranges.N);
    params.n31 = randRange(ranges.N);
    params.a2 = randRange(ranges.A);
    params.b2 = randRange(ranges.B);
    params.m2 = randRange(ranges.M);
    params.n12 = randRange(ranges.N);
    params.n22 = randRange(ranges.N);
    params.n32 = randRange(ranges.N);
end

%Show Params
function createAndDisplayShape(params, currentShape, totalShapes)
    SuperFormula3D(params.a1, params.b1, params.m1, params.n11, params.n21, params.n31, ...
                   params.a2, params.b2, params.m2, params.n12, params.n22, params.n32);
    close all;
    fprintf('Generated shape %d of %d\n', currentShape, totalShapes);
end

%Do the job
function XX = randRange(range)
    XX = range(1) + (range(2) - range(1)) * rand();
end

% Call the main function to run the program
main();
classdef KrylovPrecond

  %From Dhairya Malhotra Sept 2024

  properties

    mat_lst = {}

  end

  methods

    function p = KrylovPrecond()
    end

    function n = size(obj)
      n = length(obj.mat_lst)
    end

    function y = Apply(obj, x)
      y = x;
      n = length(obj.mat_lst)/2;
      for i = 1:n
        y = y + obj.mat_lst{(n-i)*2+1} * ( obj.mat_lst{(n-i)*2+2}' * y );
      end
    end

    function obj = Append(obj, U, Qt)
      obj.mat_lst{end+1} = U;
      obj.mat_lst{end+1} = Qt;
    end

  end

end

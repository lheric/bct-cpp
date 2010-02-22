bct_test_setup

% modularity_dir
for i = 1:size(m)(2)
	[Ci Q] = modularity_dir(m{i});
	[Ci_cpp Q_cpp] = modularity_dir_cpp(m{i});
	if length(unique(Ci)) ~= length(unique(Ci_cpp')) || unique(Ci) ~= unique(Ci_cpp')
		result = 0;
	else
		result = 1;
		map = zeros(max(Ci), 1);
		for j = 1:length(Ci)
			if ~map(Ci(i))
				map(Ci(i)) = Ci_cpp(i);
			elseif map(Ci(i)) ~= Ci_cpp(i)
				result = 0;
				break
			end
		end
	end
	bct_test(sprintf("modularity_dir %s Ci", mname{i}), result)
	bct_test(sprintf("modularity_dir %s Q", mname{i}), abs(Q - Q_cpp) < 1e-6)
end

% modularity_und
for i = 1:size(m)(2)
	ms = m{i} | m{i}';
	[Ci Q] = modularity_und(ms);
	[Ci_cpp Q_cpp] = modularity_und_cpp(ms);
	if length(unique(Ci)) ~= length(unique(Ci_cpp')) || unique(Ci) ~= unique(Ci_cpp')
		result = 0;
	else
		result = 1;
		map = zeros(max(Ci), 1);
		for j = 1:length(Ci)
			if ~map(Ci(i))
				map(Ci(i)) = Ci_cpp(i);
			elseif map(Ci(i)) ~= Ci_cpp(i)
				result = 0;
				break
			end
		end
	end
	bct_test(sprintf("modularity_und %s Ci", mname{i}), result)
	bct_test(sprintf("modularity_und %s Q", mname{i}), abs(Q - Q_cpp) < 1e-6)
end

bct_test_teardown

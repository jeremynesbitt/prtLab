% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function q = prtNorm(v)
%prtNorm normalization of vectors, whether real or complex
q = v/sqrt(sum(v.*v));
end
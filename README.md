# RbxPolyReg

Polynomial regression solver.

Full credit to the original author, Paul Lutus, and his [PolySolve](https://arachnoid.com/polysolve/) site. This is simply a port of his JavaScript code into Lua, which is covered under the GNU General Public License.

---------------

## Usage

```lua
local RbxPolyReg = require(rbxPolyRegModuleScript)

-- Get matrix functions:
local matFuncs = RbxPolyReg.MatFunctions.new()

-- Each data entry can be 1 of 3 different types:
	-- 1) {Number, Number}
	-- 2) {X = Number; Y = Number} OR {x = Number; Y = Number}
	-- 3) Vector2
local data = {
	{-1, -1},
	{0, 3},
	{1, 2.5},
	{2, 5},
	{3, 4},
	{5, 2},
	{7, 5},
	{9, 4}
}

local degrees = 2

-- Solve polynomial regression:
local result = matFuncs:ProcessData(data, degrees)
-- result[1]: Table of terms
-- result[2]: Correlation coefficient
-- result[3]: Standard error

-- List terms:
local terms = result[1]
for i,v in ipairs(terms) do
	print(("%.16e * x^%i"):format(v, i - 1))
end

-- Create source code for the polynomial regression:
local clampMin, clampMax = 0, 1000
local funcSrcCode = RbxPolyReg.ToFunctionSourceCode(terms, clampMin, clampMax)
print(funcSrcCode)

```
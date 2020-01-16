--[[
	
	Original source code is by Paul Lutus and is under the GNU General Public License:
	
		----------------------------------------------------------------
	
		Copyright (C) 2019 by Paul Lutus
		lutusp@arachnoid.com
	
		License: https://www.gnu.org/licenses/
	
		----------------------------------------------------------------
	
		See Paul Lutus' PolySolve site: https://arachnoid.com/polysolve/
	
	--------------------------------------------------------------------------------------
	
	Copyright (C) 2020 by Stephen Leitnick
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	
]]

-- This is all translated from the original JavaScript source into Lua

local RbxPolyReg = {}

------------------------------------------------------------------------------------------------------

-- The Array structure immitates 0-based index arrays to help with the port from JavaScript to Lua
local Array = {}
Array.__index = Array

function Array.new(size, default)
	return setmetatable({
		_tbl = table.create(size, default);
		__is_array_mp__ = true;
	}, Array)
end

function Array.FromTable(tbl)
	if (Array.IsArray(tbl)) then return tbl end
	local array = Array.new(#tbl)
	local newTbl = {}
	table.move(tbl, 1, #tbl, 1, newTbl)
	array._tbl = newTbl
	return array
end

function Array.IsArray(tbl)
	return (type(tbl) == "table" and rawget(tbl, "__is_array_mp__") == true)
end

function Array:Size()
	return #self._tbl
end

function Array:ForEach(func)
	for i,v in ipairs(self._tbl) do
		func(i - 1, v)
	end
end

function Array:ForEachMut(func)
	local n = #self._tbl
	for i = 1,n do
		self._tbl[i] = func(i - 1, self._tbl[i])
	end
end

function Array:__index(i)
	if (type(i) == "number") then
		return self._tbl[i + 1]
	else
		return Array[i]
	end
end

function Array:__newindex(i, v)
	self._tbl[i + 1] = v
end

------------------------------------------------------------------------------------------------------

-- A simple XY data structure
local Pair = {}
Pair.__index = Pair
RbxPolyReg.Pair = Pair

function Pair.new(x, y)
	return setmetatable({x = x or 0; y = y or 0; __is_pair_mp__ = true}, Pair)
end

function Pair.IsPair(obj)
	return (type(obj) == "table" and rawget(obj, "__is_pair_mp__") == true)
end

function Pair.Coerce(obj)
	if (typeof(obj) == "Vector2" or typeof(obj) == "Vector3") then
		-- Convert Vector2s into Pair structure:
		return Pair.new(obj.X, obj.Y)
	elseif (type(obj) == "table") then
		-- Grab table from the Array structure:
		if (Array.IsArray(obj)) then
			obj = obj._tbl
		end
		-- Copy the Pair object:
		if (Pair.IsPair(obj)) then
			return Pair.new(obj.x, obj.y)
		end
		if (#obj == 2) then
			-- {Number, Number}
			if (type(obj[1] == "number" and type(obj[2]) == "number")) then
				return Pair.new(obj[1], obj[2])
			end
		elseif (type(obj.x) == "number" and type(obj.y) == "number") then
			-- {x = Number; y = Number}
			return Pair.new(obj.x, obj.y)
		elseif (type(obj.X) == "number" and type(obj.Y) == "number") then
			-- {X = Number; Y = Number}
			return Pair.new(obj.X, obj.Y)
		end
	end
	error("Failed to coerce type \"" .. typeof(obj) .. "\" to Pair")
end

function Pair:__tostring()
	return (self.x + ", " + self.y)
end

------------------------------------------------------------------------------------------------------

-- Gauss-Jordan matrix manipulation functions
local GJ = {}
GJ.__index = GJ
RbxPolyReg.GJ = GJ

function GJ.new()
	local self = setmetatable({}, GJ)
	return self
end

function GJ:Divide(A, i, j, m)
	for q = j + 1, m - 1 do
		A[i][q] = (A[i][q] / A[i][j])
	end
	A[i][j] = 1
end

function GJ:Eliminate(A, i, j, n, m)
	for k = 0, n - 1 do
		if (k ~= i and A[k][j] ~= 0) then
			for q = j + 1, m - 1 do
				A[k][q] = (A[k][q] - (A[k][j] * A[i][q]))
			end
			A[k][j] = 0
		end
	end
end

function GJ:Echelonize(mat)
	local n = mat:Size()
	local m = mat[0]:Size()
	local i = 0
	local j = 0
	local k
	while (i < n and j < m) do
		-- Look for non-zero entries in column 'j' at or below row 'i'
		k = i
		while (k < n and mat[k][j] == 0) do
			k = (k + 1)
		end
		-- If any entry is found at row 'k'
		if (k < n) then
			-- If 'k' is not 'i', then swap row 'i' with row 'k'
			if (k ~= i) then
				mat[i], mat[k] = mat[k], mat[i]
			end
			-- If A[i][j] is not 1, divide row 'i' by A[i][j]
			if (mat[i][j] ~= 1) then
				self:Divide(mat, i, j, m)
			end
			-- Eliminate all other non-zero entries
			self:Eliminate(mat, i, j, n, m)
			i = (i + 1)
		end
		j = (j + 1)
	end
	-- Extract result column
	local terms = Array.new(m)
	mat:ForEach(function(l, v)
		terms[l] = v[n]
	end)
	return terms
end

------------------------------------------------------------------------------------------------------

local MatFunctions = {}
MatFunctions.__index = MatFunctions
RbxPolyReg.MatFunctions = MatFunctions

function MatFunctions.new()
	local self = setmetatable({
		GJ = RbxPolyReg.GJ.new();
	}, MatFunctions)
	return self
end

-- A weak substitute for printf()
function MatFunctions:NumberFormat(n, p, w)
	local s = ("%." .. tostring(p) .. "e"):format(n)
	if (#s < w) then
		s = ((" "):rep(w - #s) .. s)
	end
	return s
end

-- Produce a single 'y' result for a given 'x'
function MatFunctions:Regress(x, terms)
	local y = 0
	local exp = 1
	terms:ForEach(function(i, n)
		y = (y + (n * exp))
		exp = (exp * x)
	end)
	return y
end

-- Compute correlation coefficient
function MatFunctions:CorrelationCoefficient(data, terms)
	local n = data:Size()
	local sx, sx2, sy, sy2, sxy = 0, 0, 0, 0, 0
	local x, y
	data:ForEach(function(i, pr)
		x = self:Regress(pr.x, terms)
		y = pr.y
		sx = (sx + x)
		sy = (sy + y)
		sxy = (sxy + (x * y))
		sx2 = (sx2 + (x * x))
		sy2 = (sy2 + (y * y))
	end)
	local div = math.sqrt((sx2 - (sx * sx) / n) * (sy2 - (sy * sy) / n))
	if (div ~= 0) then
		local z = ((sxy - (sx * sy) / n) / div)
		return (z * z)
	end
	return 0
end

-- Compute standard error
function MatFunctions:StandardError(data, terms)
	local n = data:Size()
	if (n > 2) then
		local a = 0
		local q
		data:ForEach(function(_, pr)
			q = (self:Regress(pr.x, terms) - pr.y)
			a = (a + (q * q))
		end)
		return math.sqrt(a / (n - 2))
	end
	return 0
end

-- Create regression coefficients for provided data set
--    @param data: Pair array
--    @param p:    Polynomial degree
function MatFunctions:ComputeCoefficients(data, p)
	p = (p + 1)
	local n = data:Size()
	local rs = (2 * p - 1)
	-- Create square matrix with added RH column
	local m = Array.new(p)
	for i = 0, p - 1 do
		local mm = Array.new(p + 1, 0)
		m[i] = mm
	end
	-- Create array of precalculated matrix data
	local mpc = Array.new(rs, 0)
	mpc[0] = n
	data:ForEach(function(i, pr)
		-- Process precalculation array
		local x = pr.x
		local t = x
		for r = 1, rs - 1 do
			mpc[r] = (mpc[r] + x)
			x = (x * t)
		end
		-- Process RH column cells
		m[0][p] = (m[0][p] + pr.y)
		x = pr.x
		t = x
		for r = 1, p - 1 do
			m[r][p] = (m[r][p] + (x * pr.y))
			x = (x * t)
		end
	end)
	-- Populate square matrix section
	for r = 0, p - 1 do
		for c = 0, p - 1 do
			m[r][c] = mpc[r + c]
		end
	end
	-- Reduce matrix & return terms
	return self.GJ:Echelonize(m)
end

-- Test the system using known data
function MatFunctions:Test()
	local xd = {-1,0,1,2,3,5,7,9}
	local yd = {-1,3,2.5,5,4,2,5,4}
	local data = Array.new(#xd, 0)
	data:ForEachMut(function(i)
		return Pair.new(xd[i + 1], yd[i + 1])
	end)
	local terms = self:ComputeCoefficients(data, 4)
	local prec = 16
	local width = 24
	terms:ForEach(function(i, v)
		print(self:NumberFormat(v, prec, width) .. " * x^" .. i)
	end)
	local cc = self:CorrelationCoefficient(data, terms)
	local se = self:StandardError(data, terms)
	print("cc = " .. self:NumberFormat(cc, prec, width))
	print("se = " .. self:NumberFormat(se, prec, width))
end

-- Process data
--    @param data: Pair array
--    @param p:    Polynomial degree
function MatFunctions:ProcessData(data, p)
	p = math.clamp(p, 0, 18)
	data = Array.FromTable(data)
	data:ForEachMut(function(i, v) return Pair.Coerce(v) end)
	local terms = self:ComputeCoefficients(data, p)
	local cc = self:CorrelationCoefficient(data, terms)
	local se = self:StandardError(data, terms)
	return {terms._tbl, cc, se}
end

------------------------------------------------------------------------------------------------------

function RbxPolyReg.ToFunctionSourceCode(terms, clampMin, clampMax)
	terms = Array.FromTable(terms)
	local code = {"local function Solve(x)\n\treturn math.clamp(\n"}
	local size = terms:Size()
	terms:ForEach(function(i, v)
		local line = ("%.16e * math.pow(x, %i)"):format(v, i)
		if (i == (size - 1)) then
			line = (line .. ",")
		else
			line = (line .. " +")
		end
		code[#code + 1] = ("\t\t" .. line .. "\n")
	end)
	code[#code + 1] = ("\t\t%i, %i\n\t)\nend"):format(clampMin, clampMax)
	return table.concat(code)
end

------------------------------------------------------------------------------------------------------

return RbxPolyReg
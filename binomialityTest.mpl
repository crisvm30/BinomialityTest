with (StringTools):
with(LinearAlgebra):

(*
 * Networks examples that you can use to test:
 * CRN1 := [["4A","A+B"]]; #Binomial
 * CRN2 := [["S[0]+E", "ES[0]"], ["ES[0]", "S[1]+E"], ["S[1]+E", "ES[1]"], ["ES[1]", "S[2]+E"], ["S[2]+F", "FS[2]"], ["FS[2]", "S[1]+F"], ["S[1]+F", "FS[1]"], ["FS[1]", "S[0]+F"]]; # NOT Binomial
 * CRN3 := [["2A+B","C"],["C","A"]]; #Binomial
 *)


(* binomialViaHNF:
 * Procedure that receives a reversible chemical reaction network as a parameter and returns a string specifying if the network is binomial or not.
 * The parameter 'network' is a list of lists where each sub list represents a reaction. Each reaction contains 2 strings that represent the reactant and product complexes of the reaction.
 * For example, the network: 2A->B,C->D will have the next representation:[["2A","B"],["C","D"]] 
 *)
binomialViaHNF := proc(network::list)

	local m;
	m := createBinomialCoeffMatrix(network);
	print("The network is:");
	print(network);
	print("The binomial coefficient matrix of the network is:");
	print(m);
	local H:= HermiteForm(m);
	print("The Hermite Normal Form of the matrix is:");
	print(H);
	local numRows := RowDimension(H);
	local numCols := ColumnDimension(H);

	local i,j,isBinomial;
	isBinomial := true;
	for i from 1 to numRows while isBinomial do
		local numOfEntries;
		numOfEntries := 0;
		for j from 1 to numCols while isBinomial do
			if H[i][j] <> 0 then
				numOfEntries := numOfEntries + 1;
				if numOfEntries > 1 then
					# If there is more than 1 entry in a row then it is not binomial
					isBinomial := false;
				end if;
			end if;	
		end do;
	end do;

	if isBinomial then
		return "Result: The network is binomial";
	else
		return "Result: The network is NOT binomial";
	end if;
end proc:


(*createBinomialCoeffMatrix:
 * Procedure that receives a reversible chemical reaction network and returns a matrix that represents the binomial coefficient matrix  of the network.
 *)
createBinomialCoeffMatrix := proc(network::list)
	local numOfReactions := nops(network);
	local speciesSet,ODEcoeffList;
	speciesSet,ODEcoeffList := computeODECoeff(network);
	local listOfRows := [];
	
	local i,k,j;
	for i from 1 to nops(speciesSet) do
		local currentRow := [];
		#Initializing all rows with 0's
		for k from 1 to numOfReactions do
			currentRow := [op(currentRow),0]
		end do;
		#Updating the entries in the rows with the corresponding coefficient
		for j from 1 to nops(ODEcoeffList[speciesSet[i]]) do
			local currReactionNum := ODEcoeffList[speciesSet[i]][j][1];
			currentRow[currReactionNum] :=  ODEcoeffList[speciesSet[i]][j][2];
			
		end do;
		listOfRows := [op(listOfRows),currentRow];
	end do;
	
	local m:= Matrix(listOfRows);
	return m;
end proc:

(*computeODECoeff:
 * Procedure that receives a chemical reaction network and returns 2 objects: a set containing the species of the network, and a table containing the ODE coefficients for each specie 
 *)
computeODECoeff := proc(network::list)
	
	local numOfReactions := nops(network);
	local speciesSet := {};
	local ODEcoeffList := table();

	local i,j,k,y,m;
	for i from 1 to numOfReactions do
		for j from 1 to nops(network[i]) do
			
			#string with the current complex(without unnecessary blank spaces)
			local auxString := SubstituteAll(network[i][j]," ","");
			
			#list that will contain the objects of the complex
			local currSpeciesList := Split(auxString,"+"); 
			
			for k from 1 to nops(currSpeciesList) do
				y := 1;
				for m from 1 to length(currSpeciesList[k]) while IsDigit(currSpeciesList[k][m]) do
					y := m+1;					 			
				end do;
				
				#Get the specie and add it to the speciesSet
				local specie := substring(currSpeciesList[k], y .. length(currSpeciesList[k]));
				speciesSet := speciesSet union {specie};

				#Get the coefficient of the current specie in the complex
				local coef;
				coef := 1; 
				if y > 1 then
					coef := parse(substring(currSpeciesList[k], 1 .. (y-1)));
				end if;

				#If we are in the product complex then the coefficient should be negative
				if j = 2 then
					coef := -coef; 
				end if;

				#Adding the coefficient of the current specie in the ODEcoeffList
				if type(ODEcoeffList[specie],'list') then
					local notFound := true;
					local b;
					for b from 1 to nops(ODEcoeffList[specie]) while notFound do
						if ODEcoeffList[specie][b][1] = i then
							local auxList :=  ODEcoeffList[specie][b];
							 ODEcoeffList[specie] := subsop(b=NULL,hashMap[specie]);
							ODEcoeffList[specie] := [op(ODEcoeffList[specie]),[auxList[1],auxList[2]+coef] ];
							notFound := false
						end if;
					end do;
					if notFound = true then
						ODEcoeffList[specie] := [op(ODEcoeffList[specie]),[i,coef] ]
					end if;
				else
					ODEcoeffList[specie] := [ [i,coef] ]
				end if;		
			end do;
		end do;
	end do;

return speciesSet,ODEcoeffList;
end proc:



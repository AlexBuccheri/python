# README

### TODOs in This Test
Try repeating this test with parsed results - need to know if it's robust 

New test should have the structure:
* Parse result
* Parse reference
* Perform comparison
* Maybe there needs to be an exception for comparison of eigenvalues 

Set up a CMakeLists.txt to run pytest: Can do so in a single command

## TODOs in the Test Suite
* Move the parsers out of the test suite directory
* Get rid of ALL the current framework for comparing values and just replace with numpy
    * This solves our reporting issue - numpy gives way more info, for example:
* Look at replacing the top-level run routine with Ctest
  
#### Already Prototyped Here 
* Scrap the init files with a dictionary of tolerances 
* Have a generic routine for comparing values in results and reference dictionaries (I need to discuss this with python experts and see if my solution is robust enough)
* I’m not keen on the NOMAD parser for a couple of reasons - it has shit loads of dependencies and the dict keys are not intuitive.
    * Pros: 
      * Maintained by Alvin
      * No duplication of work   
      * Probably faster
    * Cons:
      * ALOT of dependencies
      * Can't disentangle from NOMAD - too much unnecessary code
      * Have to install as a submodule, then get the paths consistent (if not installing)  
      * We've already written parsers for alot of the files
      * References will break as soon as we change either a) file format (i.e. move to structured) or b) change parser
    
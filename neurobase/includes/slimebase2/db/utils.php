<?php

function _and(){	
	$args = func_get_args();
	$ret = "";
	for($i = 0; $i < sizeof($args); $i++){
		if($i > 0){
			$ret .= " AND " . encapsulateClause($args[$i]);			
		}
		else{
			$ret = encapsulateClause($args[$i]);	
		}
	}
	return $ret;
}

function encapsulateClause($clause){
	if(!(strpos($clause, " AND ") === false)){
		return "(" . $clause . ")";	
	}
	if(!(strpos($clause, " OR ") === false)){
		return "(" . $clause . ")";	
	}
	return $clause;
}

function _or($clause1, $clause2){
	$args = func_get_args();
	$ret = "";
	for($i = 0; $i < sizeof($args); $i++){
		if($i > 0){
			$ret .= " OR " . encapsulateClause($args[$i]);			
		}
		else{
			$ret = encapsulateClause($args[$i]);	
		}
	}
	return $ret;
}

function select(){
	return new Select();	
}


class Select{
	
	public $whereClause = "";
	public $limitClause = "";
	public $columnClause = "";	
	public $statement = "";
	public $table = "";
	
	function __construct(){
		$this->columnClause = "*";	
	}
	
	function columns($columns){
		$this->columnClause = "";
		for($i = 0; $i < sizeof($columns); $i++){
			$this->columnClause .= $columns[$i]. ",";
		}
		$this->columnClause = rtrim($this->columnClause, ",");
	}
	
	function execute(){
		$this->statement = " SELECT " . $this->columnClause . " FROM " . $this->table . $this->whereClause . $this->limitClause;
		return mysql_query($this->statement);
	}
	
	function innerJoin($clause){
		$this->table = 	$clause;
	}
	
	function limit($index, $num){
		$this->limitClause = " LIMIT " . $index . "," . $num;
	}
		
	function where($clause){
		$this->whereClause = " WHERE " . $clause;
	}
	
	function table(){
		$args = func_get_args();
		if(is_array($args[0])){
			$rev = array_reverse($args[0]);
			$this->table = array_pop($rev);
			$el = array_pop($rev);
			while($el != NULL){
				$this->table .= "," . $el;
				$el = array_pop($rev);
			}
		}
		else{
			$this->table = $args[0];			
		}
		
		
	}
	
}

?>
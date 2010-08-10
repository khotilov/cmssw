<? include "configure.php"?>
<html>
<head>
  <meta http-equiv="content-type" content="text/html; charset=utf-8" />
  <title>HDQM Web Inteface</title>
  <style type="text/css" title="currentStyle">
    @import "/<?=$GLOBALS['userName']?>/lib/dataTables/media/css/demo_page.css";
    @import "css/dropdown.css";
    @import "/<?=$GLOBALS['userName']?>/lib/dataTables/media/css/demo_table_jui.css";
		</style>
  <style>
    #dt_example{
      padding: 20px;
    }
  </style>
<link rel="stylesheet" href="http://ajax.googleapis.com/ajax/libs/jqueryui/1.7.2/themes/smoothness/jquery-ui.css" type="text/css" />
    <!-- @import "/<?=$GLOBALS['userName']?>/lib/dataTables/media/prettyPhoto/prettyPhoto/css/prettyPhoto.css"; -->
  </style>
  <script type="text/javascript" src="/<?=$GLOBALS['userName']?>/lib/jquery.js"></script>
  <script type="text/javascript" src="/<?=$GLOBALS['userName']?>/lib/jquery.blockUI.js"></script>

  <!-- <script type="text/javascript" src="/<?=$GLOBALS['userName']?>/lib/dataTables/media/prettyPhoto/jquery.prettyPhoto.js"></script> -->

  <!-- <script type="text/javascript" src="/<?=$GLOBALS['userName']?>/lib/dataTables/media/jquery.js"></script> -->
  <script type="text/javascript" src="/<?=$GLOBALS['userName']?>/lib/dataTables/media/js/jquery.dataTables.js"></script>
  <script type="text/javascript" src="functions.js"></script>
  <script type="text/javascript" charset="utf-8">

    var gaiSelected =  [];

	$(document).ready(function() 
	{
	
		var oTable;
		var trendsDir;
		var runTypeName="";
		var typeName="";
		
		typeName = getUrlVars()["subDet"];
		var tagVersion = getUrlVars()['subDet'];
		var tagName = "HDQM_"+typeName+"_"+tagVersion;

		
		
		$("#form").submit(function() 
		{
			if(!checkform(this))
				return false;
			$.blockUI({ message: '<h1><img src="images/bigrotation2.gif" /> Processing your request, please wait...</h1>' });

			var sData = "";
			var first = 0;
			var rows = oTable.fnGetNodes();
			for( i=0; i<rows.length; i++ ) 
			{
				var row = rows[i];
				var lastCell = row.cells.length - 1;
				if( row.cells[2].childNodes[0].checked ) 
				{
					if( first == 0 ) 
					{
						first = 1;
					}
					else 
					{
						sData += "&";
					}
					sData += $('input', oTable.fnGetNodes(i)).serialize();
					
				}
			}
			//getting the name and id for the run Type for the Plots
			var runType=document.getElementById('runtype').value;
			//concating the name 
			runTypeName=runType.substring(0,runType.lastIndexOf(':'));
			
			//getting the subDet from the URL
			typeName = getUrlVars()["subDet"];
			//getting the id for the Run Type
			var tagVersion = runType.substring(runType.lastIndexOf(':')+1,runType.length);
			var tagName = "HDQM_"+typeName+"_"+tagVersion;

			var postdata = $('#first').serialize()+"&"+$('#last').serialize()+"&runType="+runTypeName;
			alert( "The following data were submitted to the server:\n\n\nTAG NAME:"+tagName+" \n\n"+postdata+"&"+sData );

		$.ajax({
		type: "POST",
		url: '/submitForm/'+typeName+"/"+tagName,
		data: postdata+"&"+sData,
          cache: false,
		async: true,
		success: function(data) {
	    $('.result').html(data);
	    trendsDir = data.split(":")[1];
	    var len = trendsDir.length;
	    trendsDir = trendsDir.substring(2, len-2);
	   // window.open(trendsDir+'/index.html');
	  },
	  complete: function() {
	    $.unblockUI();	  
	  }
	});

        return false;
    });

    $('#subDet').html("HDQM Web Interface for " + typeName);


    //alert(typeName+'FullList.txt');
    oTable = $('#example').dataTable( {
      "iDesplayLength": 100,
      "bProcessing": true,
	  "bServerSide":false,
      //<!-- "sAjaxSource": '/dataTables/media/json_source.txt', -->
     "sAjaxSource": typeName+'FullList.txt',
	 // "sAjaxSource": "../examples_support/server_processing_id.php",
      "bStateSave": true, <!-- save the state using cookies -->
      "bJQueryUI": true,
	  "bSortClasses":false,
      "sPaginationType": "full_numbers",
      "fnRowCallback": function( nRow, aData, iDisplayIndex ) 
	  {
	  
        //additional components in each row
		//plot index
		var hiddenIndex='<input type="hidden" name="index" id="index" value="0">';        
        //log Y
		var logY = '<input type="checkbox" name="textCheck" id="logY'+iDisplayIndex+'" value="'+aData[1]+'" onclick="stateChanged()"/>';
		//selected
		var selectCheckBox = '<input type="checkbox" name="check" id="ch'+aData[1]+'" value="'+aData[1]+'" onclick="stateChanged()"/>';
    
		
        $('td:eq(1)', nRow).html(logY);
        $('td:eq(2)', nRow).html(selectCheckBox);
		$('td:eq(3)', nRow).html(hiddenIndex);
		
		
		//add the selected css class for selecting rows
		if ( jQuery.inArray(aData[1], gaiSelected) != -1 )
		{
			$(nRow).addClass('row_selected');
		}	
        return nRow;
      },
	  //columns
      "aoColumns": [
        null,
		{"sClass": "center"},
        {"sClass": "center"}

		
      ]
    });
	
	/* Click event handler */
	$('#example tbody tr').live('click', function () {
	//	alert((this));
		var aData = oTable.fnGetData( this );
		var iId = aData[1];
		
		if ( jQuery.inArray(iId, gaiSelected) == -1 )
		{			
		    if ((!stateChange)&&(this).cells[3].childNodes[0].value==0)
			{
				gaiSelected[gaiSelected.length++] = iId;
		       
			}
		}
		else
		{
			//deselect
			if (!stateChange&&(this).cells[3].childNodes[0].value==0)
			{
			gaiSelected = jQuery.grep(gaiSelected, function(value) {
				return value != iId;
			} );
			}
		}
	
	if (!stateChange&&(this).cells[3].childNodes[0].value==0)
	{
			$(this).toggleClass('row_selected');
			
	}
	else
	{
		stateChange=false;
	}
	} );
	}
	);
<!--validation check
function checkform ( form )
{

  if ((form.last.value<=form.first.value))
  {
		alert("The first run value cannot be greater than the Last Run Value");
		return false;
  }
  if (form.first.value<0) 
  {
    alert( "The first run value cannot be negative" );
    return false ;
  }
  if (form.last.value>999999999999)
  {
	alert("The last run value cannot be greater than 999999999999");
    return false ;
  }
  return true;
  }
  
  //global Variables
  var histograms=0;
  var stateChange=false;
  
    //function for checking if the selected check Box changed state
  function stateChanged()
  {
	stateChange=true;
  }

  function addToNewHistogram(index,color2)
  {	
	//histograms++;
	if (histograms>5)
	{
		alert("You can only select up to 5 different plots");
		return;
	}
	//checking if no row is selected
	if (gaiSelected.length==0)
	{
		alert("You have to select at least one row to add it to a new histogram");
		return;
	}
	
	
	
	//loop to get all the selected rows
	for (var i=0;i<gaiSelected.length;i++)
	{
		//getting the current cell
		var currentCell=document.getElementById("ch"+gaiSelected[i]);
		//alert(gaiSelected[i]);
		//alert(currentCell);
		//removing the selected class from the current row
		$(currentCell.parentNode.parentNode).removeClass("row_selected");
		//changing the color of the current Row
		currentCell.parentNode.parentNode.style.backgroundColor = color2;
		//adding the histogram index at the current row
		currentCell.parentNode.parentNode.cells[3].childNodes[0].value=index;
		//checking the checkbox
		currentCell.checked=true;
		
	
	}
	//clearing the array of the selected rows
	gaiSelected=[];

  }
 
	
	function clearAllPromt()
	{
		var answer = confirm ("Are you sure you want to clear all your selections?");
		if (answer)
			window.location.reload();
		
	}

	//Populating the RunType list from the XML file
	function loadList(type)
	{
	var x=[];
	var runTypeNamesArray=[];
	var runTypeIdsArray=[];

	var type2="SiStrip";
	
		if (window.XMLHttpRequest)
		{// code for IE7+, Firefox, Chrome, Opera, Safari
			xmlhttp=new XMLHttpRequest();
		}
		else
		{// code for IE6, IE5
			xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
		}
		xmlhttp.onreadystatechange=function()
		{
			if (xmlhttp.readyState==4 && xmlhttp.status==200)
			{
				
				//Getting the root node of the XML File
				var root=xmlhttp.responseXML.documentElement;
				var x=root.getElementsByTagName("Type");
				for (i=0;i<x.length;i++)
				{
					if ((x[i].getAttribute("name"))==(type))
					{
						xx=x[i].getElementsByTagName("RunType");
						for (j=0;j<xx.length;j++)
						{
							runTypeIdsArray[runTypeIdsArray.length++] =xx[j].getAttribute("id");
							runTypeNamesArray[runTypeNamesArray.length++]=xx[j].firstChild.nodeValue;
						}
						break;
					}//if
				}//for

			}
			var text="";			
			//creating the option list html element
			for(var i=0;i<runTypeIdsArray.length;i++)
			{
				text+="<option value= "+ runTypeNamesArray[i]+":"+runTypeIdsArray[i] + ">"+runTypeNamesArray[i]+"</option>";
			}
			$('#runtype').html(text);
		}		
		
		xmlhttp.open("GET","RunTypesData.xml",true);
		xmlhttp.send();
	}
 
</script>



</head>
<body id="dt_example" onload="loadList(getUrlVars()['subDet'])">
<div align= "left">
<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
   <a href="./" align="right">Back to Main</a><br> <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
   <a  href="javascript:clearAllPromt()">Clear All</a> 
</div>
	<h1 id="title"><div id="subDet"></div></h1>
    <p>
	
      Select the variables you want to plot and the range of runs. Once done press the submit button to create the trend plots. <br>
      - Note that you must allow pop-ups from this page to see the results -
    </p>
    <div class="demo_jui">
	<!--testing-->
	<div>
	 <ul class="dropdown">
        	<li><a href="#">Add To new Histogram</a>
        		<ul class="sub_menu">
        			<li><a href="javascript:addToNewHistogram(1,'0066CC')"><font color="0066CC">Plot 1</font></a></li>
					<li><a href="javascript:addToNewHistogram(2,'006633')"><font color="006633">Plot 2</font></a></li>
					<li><a href="javascript:addToNewHistogram(3,'CCCC33')"><font color="CCCC33">Plot 3</font></a></li>
					<li><a href="javascript:addToNewHistogram(4,'FF3333')"><font color="FF3333">Plot 4</font></a></li>
					<li><a href="javascript:addToNewHistogram(5,'6699CC')"><font color="6699CC">Plot 5</font></a></li>
				</ul>
			</li>
	</ul>
	</div>
	<!--/testing-->
	
	
      <form id="form" method="post"> 
	
<div style="text-align:center; padding-bottom:1em;">
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;PARAMETERS: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;First Run: 
	<input type="text" name="first" id="first" value="0" />
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Last Run: 
	<input type="text" name="last" id="last" value="9999999"/>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run Type:
	
<select id="runtype" name="sada">
</select>

	 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <input type="submit" value="Submit"/>
	</div>
	<table cellpadding="0" cellspacing="0" border="0" class="display" id="example">
	  <thead>
	    <tr>
	      <th width="50%">Name</th>
	      <th width="25%">LogY</th>
		  <th width="25%">Selected</th>


	    </tr>
	  </thead>
	  <tbody>
	 
	  </tbody>
	  <tfoot>
	  </tfoot>
	</table>
	<br>
      </form>
<!--    </div> -->
  </div>
<div id=processing style="display:none;">
<div style="text-align:center;padding:15px;font: normal 15px Arial,
Helvetica, sans-serif;color:#000000;font-weight:bold;width:350px">
<div class="BoxTitle" style="text-align:center;"></div>
<img src="images/bigrotation2.gif" style="margin-top:10px">
<p>Please Stand By......</p>
</div>
</div>
</body>
</html>

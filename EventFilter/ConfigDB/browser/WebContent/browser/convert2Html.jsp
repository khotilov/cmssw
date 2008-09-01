<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<%@page import="java.io.PrintWriter"%>
<%@page import="java.io.ByteArrayOutputStream"%>
<%@page import="confdb.data.IConfiguration"%>
<%@page import="confdb.converter.ConverterBase"%>
<%@page import="confdb.converter.OnlineConverter"%>
<%@page import="confdb.converter.ConverterException"%>
<%@page import="confdb.converter.BrowserConverter"%>
<%@page import="confdb.db.ConfDBSetups"%>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
<title>HLT config</title>

<script type="text/javascript" src="../js/yui/yahoo/yahoo.js"></script>
<script type="text/javascript" src="../js/yui/utilities/utilities.js"></script>
<script type="text/javascript" src="../js/yui/event/event.js"></script>

<style type="text/css">

body {
	margin:0;
	padding:0;
	border: 1px solid #B6CDE1; 
<%
  String background = request.getParameter( "bgcolor" );
  if ( background != null )
	  out.println( "background:#" + background + ";" );
%>
}

</style>


<script type="text/javascript">

function signalReady()
{
  if ( parent &&  parent.iframeReady )
    parent.iframeReady();
}

//YAHOO.util.Event.onDOMReady( signalReady );

 </script>

</head>

<body onload="signalReady()">
<pre>
<%
  try {
	String index = request.getParameter( "dbIndex" );
	if ( index == null )
	{
		String dbName = request.getParameter( "dbName" );
		if ( dbName == null ) 
		{
			out.print( "ERROR!\ndbIndex or dbName must be specified!");
			return;
		}
		else	
		{
			if ( dbName.equalsIgnoreCase( "hltdev" ) )
				dbName = "HLT Development";
			ConfDBSetups dbs = new ConfDBSetups();
		  	String[] labels = dbs.labelsAsArray();
	  		for ( int i = 0; i < dbs.setupCount(); i++ )
	  		{
	  			if ( dbName.equalsIgnoreCase( labels[i] ) )
	  			{
	  				index = "" + i;
	  				break;
	  			}
	  		}
	  		if ( index == null  )
	  		{
	  			out.print( "ERROR!\ninvalid dbName!");
	  			return;
	  		}
	  	}
	}

	int dbIndex = Integer.parseInt( index );

	ConverterBase converter = BrowserConverter.getConverter( dbIndex );

	String configName = request.getParameter( "configName" );
	String configId = request.getParameter( "configKey" );
	if ( configId == null  &&  configName == null )
	{
		out.print( "ERROR!\nconfigKey or configName must be specified!");
		return;
	}

	int configKey = ( configId != null ) ?
    	Integer.parseInt(configId) : converter.getDatabase().getConfigId(configName);

	IConfiguration conf = converter.getConfiguration( configKey );

	if ( conf == null )
		out.print( "ERROR!\nconfig " + configKey + " not found!" );
	else
	{
		String confString = null;
		try {
			if ( converter instanceof OnlineConverter )
				confString = ((OnlineConverter)converter).getEpConfigString( configKey );
			else
				confString = converter.getConverterEngine().convert( conf );
		} catch ( ConverterException e1 ) {
			BrowserConverter.clearCache();
			System.out.println( "reloading config " + configKey );
			if ( converter instanceof OnlineConverter )
				confString = ((OnlineConverter)converter).getEpConfigString( configKey );
			else
				confString = converter.getConverterEngine().convert( converter.getConfiguration( configKey ) );
		}
		out.println( confString );
	}
  } catch ( Exception e ) {
	  out.print( "ERROR!\n\n" ); 
	  ByteArrayOutputStream buffer = new ByteArrayOutputStream();
	  PrintWriter writer = new PrintWriter( buffer );
	  e.printStackTrace( writer );
	  writer.close();
	  out.println( buffer.toString() );
  }
%>
</pre>
</body>
</html>


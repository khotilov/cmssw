<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
<title>HLT config</title>

<%
  String db = request.getParameter( "db" );
  String yui = "../js/yui";
  String css = "../css";
  String js = "../js";
  String img = "../img";
  boolean online = false;
  if ( db.equals( "online" ) )
  {
	online = true;
  	yui = "../../../gui/yui";
    css = "../../css";
    js = "../../js";
    img = "../../img";
  }
%>

<link rel="stylesheet" type="text/css" href="<%=yui%>/reset-fonts-grids/reset-fonts-grids.css" />
<link rel="stylesheet" type="text/css" href="<%=yui%>/container/assets/skins/sam/container.css" />
<link rel="stylesheet" type="text/css" href="<%=yui%>/resize/assets/skins/sam/resize.css"" />
<link rel="stylesheet" type="text/css" href="<%=yui%>/treeview/assets/skins/sam/treeview.css" />
<link rel="stylesheet" type="text/css" href="<%=yui%>/tabview/assets/tabview-core.css" />
<link rel="stylesheet" type="text/css" href="<%=yui%>/button/assets/skins/sam/button.css">
<link rel="stylesheet" type="text/css" href="<%=css%>/folders/tree.css">
<link rel="stylesheet" type="text/css" href="<%=css%>/confdb.css" />

<script type="text/javascript" src="<%=yui%>/utilities/utilities.js"></script>
<script type="text/javascript" src="<%=yui%>/cookie/cookie-min.js"></script>
<script type="text/javascript" src="<%=yui%>/datasource/datasource-beta-min.js"></script>
<script type="text/javascript" src="<%=yui%>/resize/resize-min.js"></script>
<script type="text/javascript" src="<%=yui%>/json/json-min.js"></script>
<script type="text/javascript" src="<%=yui%>/treeview/treeview-min.js"></script>
<script type="text/javascript" src="<%=yui%>/tabview/tabview.js"></script>
<script type="text/javascript" src="<%=yui%>/container/container-min.js"></script>
<script type="text/javascript" src="<%=yui%>/button/button-min.js"></script>
<script type="text/javascript" src="<%=js%>/HLT.js"></script>


<style>

body, #doc3, #pg, .blindTable { 
    padding:0px; 
    margin:0px; 
}

body {
    overflow: hidden;
    position: fixed;
}

.blindTable td {
  border:0px;
}


#mainLeft { 
	margin:0px; 
	padding: 2px 5px 0px 1px;
}

#mainRight { 
	margin:0px; 
	padding: 2px 1px 0px 0px;
}


#leftHeaderDiv { 
	margin:0px; 
	padding:0.4em; 
	border-bottom: 0px;
	height: 1.2em;
}

#rightHeaderDiv {
	border-width: 1px;
	border-style: solid;
	border-bottom: 0px;
	padding:2px; 
}


.yui-nav,
.yui-content {
	border-width: 1px; 
	border-style: solid;
}

.yui-nav {
	border-top: 0px; 
	border-bottom: 0px; 
}

#treeDiv {
    overflow: auto;
}

</style>

<script type="text/javascript">

<%
  String height = request.getParameter( "height" );
  if ( height == null )
	  out.println( "var displayHeight = 0;" );
  else
	  out.println( "var displayHeight = " + height + ";" );
  String width = request.getParameter( "width" );
  if ( width == null )
	  out.println( "var displayWidth = 0;" );
  else
	  out.println( "var displayWidth = " + width + ";" );

  out.println( "var imgDir = '" + img + "';" );
  
  if ( db != null )
  {
	  if ( online )
	  {
		out.println( "var dbName = 'online';" );
		String tabId = request.getParameter( "tabId" );
		if ( tabId == null )
		      out.println( "var pageId = 'online';" );
		else
		      out.println( "var pageId = '" + tabId + "';" );
		out.println( "var onlineMode = true;" );
		out.println( "AjaxInfo._path = '../../jsp/hlt/AjaxInfo.jsp';" );
	  }
	  else
	  {
		out.println( "var dbName = '" + db + "';" );
		out.println( "var onlineMode = false;" );
		out.println( "var pageId = '" + db + "';" );
		out.println( "AjaxInfo._path = '../browser/AjaxInfo.jsp';" );
	  }
  }
  else
  {
	  out.println( "var onlineMode = false;" );
	  out.println( "var pageId = \"any\";" );
	  out.println( "var dbName = null;" );
	  out.println( "AjaxInfo._path = '../browser/AjaxInfo.jsp';" );
  }
%>


var configFrameUrl,
	configKey,
	fullName,
	Dom = YAHOO.util.Dom,
    Event = YAHOO.util.Event,
    mainLeft = null,
    mainRight = null,
    activeMainDiv,
    displayWidth,
    resize,
    oldWidth = "200px",
    detailsMode = true,
    hltCookie = null,
    cookieExpires,
    tree,
  	tooltip, 
  	tooltipElements = [],
  	tabView, 
  	tabReady = [],
    activeTab = 1;
var filter = '/cdaq/.*';
	
function init() 
{
    cookieExpires = new Date();
    cookieExpires.setFullYear( cookieExpires.getFullYear() + 1 );
    if ( parent && parent.cookie )
    {
      hltCookie = YAHOO.util.Cookie.getSubs( pageId );
      if ( hltCookie == null )
    	hltCookie = new Object();
    }

    mainLeft = Dom.get('mainLeft');
    mainRight = Dom.get('mainRight');

	if ( displayHeight == 0 )
		displayHeight = Dom.getViewportHeight();
	if ( displayWidth == 0 )
		displayWidth = Dom.getViewportWidth();

    Dom.setStyle( 'collapseDiv', 'visibility', 'hidden' );
    Dom.setStyle( 'expandDiv', 'visibility', 'hidden' );
    Dom.setStyle( mainRight, 'visibility', 'hidden' );
    Dom.setStyle(  'doc3', 'height',  displayHeight + 'px' );
    Dom.setStyle(  'pg', 'height',  displayHeight + 'px' );
    Dom.setStyle(  'pg', 'width',  displayWidth + 'px' );

	var treeHeight = displayHeight - 30;
    Dom.setStyle( 'treeDiv', 'max-height',  treeHeight + 'px' );

    resize = new YAHOO.util.Resize('mainLeft', {
            proxy: true,
            ghost: true,
            handles: ['r'],
            maxWidth: displayWidth,
            minWidth: 10
        });
    resize.on('startResize', function(ev) {
//		if ( YAHOO.env.ua.ie > 0 ) 
//  		  Dom.setStyle( 'treeDiv', 'visibility', 'hidden' );
        });

    resize.on('resize', function(ev) {
            var w = ev.width;
            if ( hltCookie && w > 10 )
            {
              hltCookie.treeWidth = w;
		  	  YAHOO.util.Cookie.setSubs( pageId, hltCookie, { expires: cookieExpires } );
		  	}
            var width = displayWidth - w - 8;
            Dom.setStyle( mainRight, 'height', displayHeight + 'px' );
            Dom.setStyle( mainRight, 'width', width + 'px');
            Dom.setStyle( 'rightHeaderBottomDiv', 'width', (width - 200) + 'px' );
    });

  var treeWidth = displayWidth / 3;
  if ( hltCookie && hltCookie.treeWidth )
	  	treeWidth = hltCookie.treeWidth;
  resize.resize(null, displayHeight, treeWidth, 0, 0, true);


  //handler for expanding all nodes
  Event.on("expand", "click", function(e) {
			tree.expandAll();
			Event.preventDefault(e);
		});
		
  //handler for collapsing all nodes
  Event.on("collapse", "click", function(e) {
			tree.collapseAll();
			Event.preventDefault(e);
		});

  //handler for collapseDiv
  Event.on("collapseDiv", "click", function(e) {
			oldWidth = Dom.getStyle( mainLeft, 'width' );
	        resize.resize( null, displayHeight, 1, 0, 0, true);
            Dom.setStyle( mainLeft, 'visibility', 'hidden' );
            Dom.setStyle( 'expandDiv', 'visibility', 'visible' );
		});

  //handler for expandDiv
  Event.on( "expandDiv", "click", function(e) {
            Dom.setStyle( 'expandDiv', 'visibility', 'hidden' );
	        resize.resize( null, displayHeight, oldWidth, 0, 0, true);
            Dom.setStyle( mainLeft, 'visibility', 'visible' );
		});

  if ( hltCookie && hltCookie.activeTab )
    activeTab = hltCookie.activeTab;
  if ( !activeTab )
    activeTab = 1;
  tabView = new YAHOO.widget.TabView( 'tabView', { activeIndex : activeTab } );
  tabView.set( 'activeIndex', activeTab );

  tabView.on( 'activeTabChange', function( eventInfo ) {
	  activeTab = tabView.get( 'activeIndex' );
	  if ( !tabReady[ activeTab ] )
	    loadTab();
      if ( hltCookie )
      {
        hltCookie.activeTab = activeTab;
	    YAHOO.util.Cookie.setSubs( pageId, hltCookie, { expires: cookieExpires } );
	  }
	} );


  if ( onlineMode )
  {
    var submitButton = new YAHOO.widget.Button( { label: "Submit", id: "submitbutton", container: "buttonTD" });
    submitButton.on("click", onSubmitClick ); 	
  
    Dom.setStyle( 'buttonTD', 'visibility', 'hidden' );
    Dom.setStyle( 'downloadTD', 'visibility', 'hidden' );
  }

  if ( pageId == 'online' )
	  filter = '';
  if ( hltCookie && hltCookie.filter )
	  	filter = hltCookie.filter;
  AjaxInfo.getTree( dbName, filter, createTree );	
}
	

function onSubmitClick( event ) 
{ 
  if ( parent && parent.submitConfig )
    parent.submitConfig( configKey, fullName );
} 
	

	
function createTree( treeData )
{
	if ( treeData.exceptionThrown )
	{
		alert( treeData.exception + ': ' + treeData.message );
		return;
	}

    if ( treeData.ajaxFailure )
    {
	  alert( 'Ajax failure: ' + treeData.ajaxFailure );
	  return;
    }

	tree = new YAHOO.widget.TreeView("treeDiv");
	var parentNode = tree.getRoot();
	createTreeRecursiveLoop( parentNode, treeData );
	tree.render();
	tree.subscribe( "clickEvent", configSelected );
	var header = '<table><tr>';
	if ( filter != '' )
		header += "<td><b>filter: "  + filter + "</b></td><td><div style='width:50px'></div></td>";
	header += '<td><a id="expand" href="#">Expand all</a> <a id="collapse" href="#">Collapse all</a></td></tr></table>';
  	Dom.get( 'leftHeaderDiv' ).innerHTML = header; 
  	Dom.setStyle( 'collapseDiv', 'visibility', 'visible' );
  	
  	// uses too much CPU power!
	//tooltip = new YAHOO.widget.Tooltip( "tooltip", { context: tooltipElements } ); 

	if ( hltCookie )
  	{
    	var config = hltCookie.selectedConfig;
    	if ( config != null )
    	{ 
      	  var node = tree.getRoot();
      	  var subdirs = config.split( "/" );
      	  if ( subdirs.length > 1 && subdirs[0] == "" )
      	  {
        	subdirs.shift();
        	subdirs[0] = '/' + subdirs[0];
      		findNode( node, config, subdirs );
      	  }
    	}
  	}
}
	

function createTreeRecursiveLoop( parentNode, treeData )
{
	for ( var i = 0; i < treeData.configs.length; i++ )
	{    
		var config = treeData.configs[i];
		config.nodeData.expanded = false;
		var configNode = new YAHOO.widget.ConfigNode( config.nodeData, parentNode );
		if ( config.nodeData.title )
		  tooltipElements.push( configNode.labelElId );
		for ( var ii = 0; ii < config.subnodes.length; ii++ )
		{    
		  config.subnodes[ii].nodeData.expanded = false;
		  var subnode = new YAHOO.widget.ConfigNode( config.subnodes[ii].nodeData, configNode );
		  if ( config.subnodes[ii].nodeData.title )
		    tooltipElements.push( subnode.labelElId );
		}
	}

	for ( var i = 0; i < treeData.dirs.length; i++ )
	{    
		var dir = treeData.dirs[i];
	    var name = dir.name;
		var dirNode = new YAHOO.widget.TextNode( { label: name, expanded: false }, parentNode );
		createTreeRecursiveLoop( dirNode, dir );
	}
}



function findNode( node, config, subdirs )
{  
  if ( !node.hasChildren() )
    return;
  if ( subdirs.length == 0 )
    return;
  var nodes = node.children;
  for ( var i = 0;  i < nodes.length; i++ )
  {
    var fullName = nodes[i].data.fullName;
    if ( fullName && fullName == config )
      return;
  }

  var subdir = subdirs.shift();
  for ( var i = 0;  i < nodes.length; i++ )
  {
    var label = nodes[i].label;
    if ( label == subdir )
    {
      nodes[i].expand();
      findNode( nodes[i], config, subdirs );
    }
  }
}	
	
	
function configSelected( event )
{
  var node = event.node;
  if ( !node.data.key )
  	return;
  	
  node.focus();
  
//  if ( parent &&  parent.configSelected )
//    parent.configSelected( node.data );

  Dom.setStyle( mainRight, 'visibility', 'visible' );

  configKey = node.data.key;
  fullName = node.data.fullName;
  Dom.get( 'fullNameTD' ).innerHTML = "<b>" + fullName + "</b>";
  var fileName = node.data.name.replace( '//s/g', '_' ) + "_V" + node.data.version;

  if ( !onlineMode )
  {
    Dom.get( 'downloadTD' ).innerHTML = 'download ' 
      + '<a href="' + fileName + '.cfg?configId='+ configKey + '&dbName=' + dbName + '">cfg</a> '
      + '<a href="' + fileName + '.py?format=python&configId='+ configKey + '&dbName=' + dbName + '">py</a>';
  }
  else
  {
    Dom.setStyle( 'buttonTD', 'visibility', 'visible' );
  }

  tabReady = [];
  loadTab();
  if ( hltCookie )
  {
    hltCookie.selectedConfig = fullName;
    YAHOO.util.Cookie.setSubs( pageId, hltCookie, { expires: cookieExpires } );
  }
  return false;
}

function loadTab()
{  
  var tabDiv = 'tab' + activeTab + 'Div';
  activeMainDiv = tabDiv + 'Main';
  Dom.setStyle( activeMainDiv, 'visibility', 'hidden' );
  Dom.get( 'rightHeaderBottomDiv' ).innerHTML = '<img src="' + imgDir + '/wait.gif">';
  var xy = Dom.getXY( tabDiv );
  var height = displayHeight - xy[1] - 2;
  Dom.setStyle( activeMainDiv, 'height', height + 'px' );
  Dom.setStyle( tabDiv, 'height', height + 'px' );

  var tabContent = '<iframe src="' + buildIFrameUrl() + '" name="configIFrame" id="configFrame" width="100%" height="'+ height + '" frameborder="0" ' + (detailsMode ? '' : ' scrolling="no"') + '></iframe>';
  Dom.get( activeMainDiv ).innerHTML = tabContent;
  tabReady[ activeTab ] = true;
}

function buildIFrameUrl()
{  
  if ( activeTab == 1 )
    return "convert2Html.jsp?configKey=" + configKey + "&dbName=" + dbName + (onlineMode ? "&online=true" : ""); 
  if ( activeTab == 2 )
    return "showSummary.jsp?configKey=" + configKey + "&dbName=" + dbName + (onlineMode ? "&online=true" : "");
  return "";
}     

  
  
function iframeReady()
{
  Dom.get( 'rightHeaderBottomDiv' ).innerHTML = "";
  Dom.setStyle( activeMainDiv, 'visibility', 'visible' );
}
	
	
YAHOO.widget.ConfigNode = function(oData, oParent ) 
{
  //this.labelStyle = "icon-gen";
  //this.href = "javascript:dummy()";
  //oData.href = "javascript:dummy()";
  YAHOO.widget.ConfigNode.superclass.constructor.call( this, oData, oParent );
};

YAHOO.extend(YAHOO.widget.ConfigNode, YAHOO.widget.TextNode, 
{
  configNode: true,

  updateIcon: function() {
        if (this.hasIcon) {
            var el = this.getToggleEl();
            if (el) {
                el.className = this.getStyle();
            }
        }
    },
    
   getDepthStyle: function(depth) {
     if ( !this.hasChildren(false) && depth >= this.depth - 1 )
       return "ygtvblankdepthcell";
     else
       return (this.getAncestor(depth).nextSibling) ? 
            "ygtvdepthcell" : "ygtvblankdepthcell";
    },

    /**
     * Returns the css style name for the toggle
     * @method getStyle
     * @return {string} the css class for this node's toggle
     */
  getStyle: function() 
  {
     // location top or bottom, middle nodes also get the top style
     var loc = (this.nextSibling) ? "t" : "l";

     // type p=plus(expand), m=minus(collapase), n=none(no children)
     var type = "n";
     if ( this.hasChildren(false) )
     {
       if ( this.expanded )
         return "xygtv" + loc + "m";
       else
         return "xygtv" + loc + "p";
     }
     return "ygtv" + loc + type;
  }
});


function dummy( node )
{
}

//When the DOM is done loading, we can initialize our TreeView
//instance:
YAHOO.util.Event.onContentReady( "doc3", init );
	
</script>

</head>
<body class="yui-skin-sam">


<div id="doc3">
  <div id="pg" class="skin1">
    <div class="yui-g" id="pg-yui-g">
	  <div class="yui-u first" id="mainLeft">
    	<div id="leftHeaderDiv" class="tree1"><img src="<%=img%>/wait.gif"></div>
        <div style="position:absolute; right:8px; top:3px; z-index:1; cursor:pointer" id="collapseDiv" ><img src="<%=img%>/collapse.gif"></div>
        <div style="position:absolute; left:0px; top:2px; z-index:2; cursor:pointer;" id="expandDiv" ><img src="<%=img%>/tree/expand.gif"></div>
        <div align="left" id="treeDiv" class="tree1" style="border-top:0px;"></div>
   	  </div>

      <div class="yui-u skin1" id="mainRight">
        <div id="rightHeaderDiv" class="header1">
  		  <table width="100%" class='blindTable'>
  		  <tr>
  		    <td><table class='blindTable'><tr>
  		     <td id='fullNameTD'><b>/PATH/CONFIG/VERSION</b></td>
  		     <td align="right" id='buttonTD' style="padding-left:40px"></td>
  		    </tr></table></td>
  			<td align="right" id='downloadTD'>download</td>
		  </tr>
		  </table>
        </div>
		<div id="tabView" class="yui-navset header1">
		  <ul class="yui-nav" id="tabViewHeader">
		    <li class="disabled"><div id="rightHeaderBottomDiv"></div><a href="#tab0Div"></a></li>
		    <li><a href="#tab1Div"><em>details</em></a></li>
		    <li><a href="#tab2Div"><em>summary</em></a></li>
		  </ul>            
		  <div class="yui-content">
			 <div id="tab0Div">tab not loaded</div>
			 <div id="tab1Div" class="tab1">
			   <div id="tab1DivMain"></div>
			 </div>
			 <div id="tab2Div" class="tab1">
			   <div id="tab2DivMain"></div>
			 </div>
		  </div>
		</div>
      </div>
    </div>
  </div>
</div>
</body>
</html>


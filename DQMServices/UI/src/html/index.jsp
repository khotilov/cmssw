<%@ page import="java.sql.*,java.util.*,java.util.regex.*,cms.dqm.workflow.*,java.security.cert.*" %>
<%@ taglib prefix="dqm" uri="/WEB-INF/cmsdqmworkflow.tld" %>

<%
  
  try {
    DBWorker db = new DBWorker();
    db.executeSQL("select 1 from dual");
  } catch (Exception e) {
    out.println("No connection to the database.");
  }

  String logged_in = WebUtils.getLoggedUser(request);
  String logged_roles = "";


  Vector<String> v = WebUtils.getLoggedRoles(request);
  for (Iterator i = v.iterator(); i.hasNext(); ) {
    logged_roles = logged_roles + i.next();
    if (i.hasNext()) {
      logged_roles = logged_roles + ", ";
    }
  }

  if (logged_roles.equals("")) logged_roles = "NONE";

%>

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>

  <meta http-equiv="content-type" content="text/html; charset=utf-8" />
  <title>CMS DQM Run Registry</title>

  <link rel="stylesheet" type="text/css" href="media/flexigrid/css/flexigrid/flexigrid.css" />
  <link rel="stylesheet" type="text/css" href="media/jquery.autocomplete.css" />
  <link rel="stylesheet" type="text/css" href="media/index.css" />
  <link rel="stylesheet" type="text/css" href="media/flora.tabs.css" />

  <script type="text/javascript" src="media/jquery-1.2.6.js"></script>
  <script type="text/javascript" src="media/jquery.cookie.js"></script>
  <script type="text/javascript" src="media/jquery.timers.js"></script>
  <script type="text/javascript" src="media/jquery.menu.js"></script>
  <script type="text/javascript" src="media/jquery.tooltip.js"></script>
  <script type="text/javascript" src="media/jquery.autocomplete.js"></script>
  <script type="text/javascript" src="media/jquery.dimensions.js"></script>
  <script type="text/javascript" src="media/flexigrid/flexigrid.js"></script>
  <script type="text/javascript" src="media/utils.js"></script>
  <script type="text/javascript" src="media/ui.core.js"></script>
  <script type="text/javascript" src="media/ui.tabs.js"></script>
  <script type="text/javascript" src="media/ui.progressbar.js"></script>

<script type="text/javascript">

  function checkSession() {
    $.ajax({
      type: "GET",
      url: "heartbeat",
      dataType: "text",
      success: function(ret) {
        if (ret == '0') {
          logoutUser();
        }
      }
    });
  }

  function logoutUser() {
    frames['logout'].location.href = "https://login.cern.ch/adfs/ls/?wa=wsignout1.0";
    setTimeout("location.href=location.href.replace('https:','http:');", 200);
  }

  var total = 0;

  function getTag() {
    var tag = '';
    do { 
      tag = prompt('Please type in Tag value for COMPLETED run(s) (required)');
    } while (tag == '' || tag == null);
    return tag;
  }

  function toggleRows(number) {
    
    var plus = $("#row" + number + " TD[field=RUN_COUNT_TAGS] div").find("img[src*=plus]");
    var minus = $("#row" + number + " TD[field=RUN_COUNT_TAGS] div").find("img[src*=minus]");

    if ($("[tagParent=" + number + "]").length == 0) {
      
      var load = $("#row" + number + " TD[field=RUN_COUNT_TAGS] div").find("img[src*=tagloader]").show();
      $(plus).hide();
      $(load).show();
      
      $.ajax({
        type: "GET",
        url: "runversiondata?qtype=RUN_NUMBER&query=" + number,
        processData: false,
        dataType: "xml",
        error: function(o) {
          $(load).hide();
          $(plus).show();
          errorMessage(o);
        },
        success: function(ret) {
          var xml = $(ret).find("RESULTS");
					var rownum = 0;
          $(xml).children().each(function (i, r) {
            var row = $("#row" + number).clone();
            $(row).removeAttr("id");
            $(row).addClass("prevTag");
            $(row).removeClass("erow");
            $(row).attr("tagParent", number);
            $(r).children().each(function (i, o) {
              var t = $(o).text();
              if (o.nodeName == "RUN_NUMBER") t = "";
              else if (o.nodeName == "RUN_EVENTS" || o.nodeName == "RUN_RATE" || o.nodeName == "RUN_BFIELD") t = formatNumber(t);
              $(row).find("TD[field = " + o.nodeName + "] div").text(t);
            });

            $.each([<dqm:listSubsystemsJS/>], function(si, so) {
              var v = $(r).find(so).text();
              if (v) {
                var c = $(r).find(so + "_comment").text();
                var arrows = "&nbsp;";
                var statTip = "";
                if (c != "") {
                  arrows = "!";
                  statTip = " statTip";
                }
                var t = "<span class=\"stat stat" + v + statTip + "\" status=\"" + v + "\"\>" + arrows + "</span>";
                $(row).find("TD[field = " + so + "] div").html(t)

                if (c != "") {
                  $(row).find("TD[field=" + so + "] div").tooltip({
                    delay: 0,
                    fade: 250,
                    track: true,
                    bodyHandler: function() {
                      var comment = c;
                      return comment;
                    }
                  });
                }

              }
            });

            $(row).find("TD[field=NR] div").html("");
						if (rownum == 0) {
							$(row).find("TD[field=RUN_COUNT_TAGS]").addClass("tree_end");
						} else {
							$(row).find("TD[field=RUN_COUNT_TAGS]").addClass("tree_middle");
						}
						rownum++;
            $(row).insertAfter("#row" + number);

          });
          $(load).hide();
          $(minus).show();
        }
      });
    } else {
      $("[tagParent=" + number + "]").toggle();
      $(plus).toggle();
      $(minus).toggle();
    }

  };

  $(document).ready( function () {


    var subsystems = [<dqm:listSubsystemsJS/>];

    var preProcess = function (data) {

      if (!data) return data;

      var st = (data.page - 1) * data.rp;
      total = data.total;

      $.each(data.rows, function(i, row) {

        var number = parseInt(row["RUN_NUMBER"]);

        row.id = number;
        row["RUN_EVENTS"] = formatNumber(row["RUN_EVENTS"]);
        row["RUN_RATE"] = formatNumber(row["RUN_RATE"]);
        row["RUN_BFIELD"] = formatNumber(row["RUN_BFIELD"]);
        row["NR"] = data.total - (st + i);

        var tags = parseInt(row["RUN_COUNT_TAGS"]) - 1;
        if (tags > 0) {
          row["RUN_COUNT_TAGS"] = "<img onclick=\"toggleRows(" + number + ");\" src=\"media/plus.png\"/><img onclick=\"toggleRows(" + number + ");\" style=\"display: none\" src=\"media/minus.png\"/><img style=\"display: none\" src=\"media/tagloader.gif\"/>";
        } else {
          row["RUN_COUNT_TAGS"] = "";
        }

        $.each(subsystems, function(si, so) {
          if ((row[so]) && row[so] != "null") {
            var arrows = "&nbsp;";
            var statTip = "";
            if ((row[so + "_comment"]) && (row[so + "_comment"] != "null")) {
              arrows = "!";
              statTip = " statTip";
            }
            row[so] = "<span class=\"stat stat" + 
              row[so] + statTip + 
              "\" status=\"" + row[so] + 
              "\" comment=\"" + row[so + "_comment"] + "\">" + 
              arrows + "</span>";
          }
        });

      });

      return data;
    };

    var postProcess = function () {
      $("span.statTip").parent().tooltip({
        delay: 0,
        fade: 250,
        track: true,
        bodyHandler: function() {
          var comment = $("span.statTip", this).attr("comment");
          return comment;
        }
      });
    };

    var timerToggle = function () {
      $("#flex1").stopTime();
      if ($.cookie("flex_timer") == "true") {
        $(this).find("span").text("Manual Refresh");
        $.cookie("flex_timer", null);
      } else  {
        $(this).find("span").text("Auto Refresh");
        $("#flex1").everyTime(60 * 1000, function(i) {
          $(this).flexReload();
        });
        $.cookie("flex_timer", "true");
      }
    };

    var selectToggle = function () {
      if ($.cookie("flex_multiselect") == "true") {
        $(this).find("span").text("Single Select");
        $("#flex1").flexOptions({ singleSelect: true });
        $.cookie("flex_multiselect", null);
      } else  {
        $(this).find("span").text("Multi Select");
        $("#flex1").flexOptions({ singleSelect: false });
        $.cookie("flex_multiselect", "true");
      }
    };

    var tableToggle = function () {
      if ($.cookie("flex_table_summary") == "true") {
        $("#flex1").flexOptions({ url: "runregisterdata", newp: 1 });
        $("#flex1").flexReload();
        $(this).find("span").text("Run Registry Data");
        $.cookie("flex_table_summary", null);
      } else {
        $("#flex1").flexOptions({ url: "runinfodata", newp: 1 });
        $("#flex1").flexReload();
        $(this).find("span").text("All Run Info Data");
        $.cookie("flex_table_summary", "true");
      }
    }

    var canEdit = function (status) {
      if (status == "COMPLETED") return false;
      if (<%= WebUtils.hasLoggedRole(request, WebUtils.EXPERT) %>) return true;
      if (status == "ONLINE" && <%= WebUtils.hasLoggedRole(request, WebUtils.ONLINE) %>) return true;
      if (status == "OFFLINE" && <%= WebUtils.hasLoggedRole(request, WebUtils.OFFLINE) %>) return true;
      return false;
    }

    var editPress = function () {

      var number = parseInt($("div.button_edit").attr("run_number"));

      if (!number) {
        alert("Run not selected.");
        return;
      }

      $.showRunEditForm(number);
    }

    var goLink = function (l) {
      var number = parseInt($("div.button_edit").attr("run_number"));

      if (!number) {
        alert("Run not selected.");
        return;
      }

      l = l.replace(/\{RUN_NUMBER\}/g, number);
      window.open(l);
    }

    var onRowSelected = function (row) {

      var number = parseInt($("td[field=RUN_NUMBER] div", row).text());
      var status = $("td[abbr=RUN_STATUS] div", row).text();

      if (trim(status).length < 2) status = "ONLINE";

      $("div.button_edit").attr("run_number", number);

      var allow = canEdit(status);
      if (allow) {
        if ($("td[abbr=RUN_EVENTS] div", row).text() == "null" || 
            $("td[abbr=RUN_START_TIME] div", row).text() == "" ||
            $("td[abbr=RUN_RATE] div", row).text() == "null" ||
            $("td[abbr=RUN_L1KEY] div", row).text() == "" ||
            $("td[abbr=RUN_HLTKEY] div", row).text() == "") {
          allow = 0;
        }
      }

      if (allow) {
        $("div.button_edit span").text("Edit #"+ number);
      } else {
        $("div.button_edit span").text("View #"+ number);
      }

    }

    var onRowDblClick = function (row) {
      var number = parseInt($("td[field=RUN_NUMBER] div", row).text());
      $("div.button_edit").attr("run_number", number);
      editPress();
    }

    /* DataTable default preferences */
    var fp = {
      url: 'runregisterdata?',
      dataType: 'json',
      method: 'get',
      colModel : [
        {display: 'ID', name : 'RUN_ID', width : 20, sortable : true, align: 'right', hide: 'true'},
        {display: 'Tags', name : 'RUN_COUNT_TAGS', width : 20, sortable : true, align: 'center'},
        {display: '#', name : 'NR', width : 40, align: 'right', hide: 'true'},
        {display: 'Run#', name : 'RUN_NUMBER', width : 35, sortable : true, align: 'right'},
        {display: 'Global Name', name : 'RUN_GLOBALNAME', width : 80, sortable : true, align: 'left', hide: 'true'},
        {display: 'Online Shifter', name : 'RUN_ONLINE_SHIFTER', width : 80, sortable : true, align: 'left'},
        {display: 'Offline Shifter', name : 'RUN_OFFLINE_SHIFTER', width : 80, sortable : true, align: 'left', hide: 'true'},
        {display: 'Status', name : 'RUN_STATUS', width : 100, sortable : true, align: 'center'},
        {display: 'Current Tag', name : 'RUN_TAG', width : 70, sortable : true, align: 'left'},
        {display: 'Started', name : 'RUN_START_TIME', width : 100, sortable : true, align: 'left'},
        {display: 'Stopped', name : 'RUN_END_TIME', width : 100, sortable : true, align: 'left', hide: 'true'},
        {display: 'Events', name : 'RUN_EVENTS', width : 50, sortable : true, align: 'right'},
        {display: 'Rate, hz', name : 'RUN_RATE', width : 50, sortable : true, align: 'right'},
        {display: 'B Field', name : 'RUN_BFIELD', width : 50, sortable : true, align: 'right'},
        {display: 'B Field Comment', name : 'RUN_BFIELD_COMMENT', width : 130, sortable : true, align: 'left', hide: 'true'},
        {display: 'Stop Reason', name: 'RUN_STOP_REASON', width : 130, sortable : true, align: 'left'},
        {display: 'L1 Key', name: 'RUN_L1KEY', width : 130, sortable: true, align: 'left'},
        {display: 'HLT Key', name: 'RUN_HLTKEY', width : 130, sortable: true, align: 'left'},
        {display: 'L1 Sources', name: 'RUN_L1SOURCES', width : 150, sortable: true, align: 'left', hide: 'true'},
        {display: 'Online Comment', name: 'RUN_ONLINE_COMMENT', width : 150, sortable: true, align: 'left'},
        {display: 'Offline Comment', name: 'RUN_OFFLINE_COMMENT', width : 150, sortable: true, align: 'left', hide: 'true'},
        <dqm:listSubsystems type="ONLINE">
        {display: '${sub_abbr}', name: '${sub_abbr}', width: 30, sortable: true, align: 'center'},
        </dqm:listSubsystems>
        <dqm:listSubsystems type="OFFLINE">
        {display: '${sub_abbr}', name: '${sub_abbr}', width: 30, sortable: true, align: 'center', hide: 'true'},
        </dqm:listSubsystems>
        {display: 'Created', name: 'RUN_CREATE_TIME', width: 100, sortable: true, align: 'left', hide: 'true'}
      ],
      preProcess: preProcess,
      postProcess: postProcess,
      onRowSelected: onRowSelected,
      onRowDblClick: onRowDblClick,
      autoload: false,
      searchitems : [
        {display: 'RUN#', name : 'RUN_NUMBER', isdefault: true},
        {display: 'Online Shifter', name : 'RUN_ONLINE_SHIFTER'},
        {display: 'Offline Shifter', name : 'RUN_OFFLINE_SHIFTER'},
        {display: 'Global Name', name : 'RUN_GLOBALNAME'},
        {display: 'Status', name : 'RUN_STATUS'},
        {display: 'Current Tag', name : 'RUN_TAG'},
        {display: 'Tags', name : 'RUN_COUNT_TAGS'},
        {display: 'L1 Key', name: 'RUN_L1KEY'},
        {display: 'HLT Key', name: 'RUN_HLTKEY'},
        {display: 'Started', name : 'RUN_START_TIME'},
        {display: 'Stopped', name : 'RUN_END_TIME'},
        {display: 'Events', name : 'RUN_EVENTS'},
        {display: 'Rate, hz', name : 'RUN_RATE'},
        <dqm:listSubsystems type="ONLINE">
        {display: '${sub_abbr}', name: '${sub_abbr}'},
        </dqm:listSubsystems>
        <dqm:listSubsystems type="OFFLINE">
        {display: '${sub_abbr}', name: '${sub_abbr}'},
        </dqm:listSubsystems>
        {display: 'Created', name : 'RUN_CREATE_TIME'}
      ],
      sortname: "RUN_NUMBER",
      sortorder: "desc",
      rp: 20,
      rpOptions: [10,15,20,25,30,35,40,45,50],
      height: 480,
      savePrefix: 'flex',
      singleSelect: true,
      width: 'auto',
      usepager: true,
      useRp: true,
      showTableToggleBtn: true,
      title: '',
      errormsg: 'Database error',
      buttons : [
        {name: 'Manual Refresh', bclass: 'refresh', dclass: 'button_refresh', onpress : timerToggle },
        {name: 'Single Select', dclass: 'button_select', onpress : selectToggle },
        {separator: true},
        {name: 'Columns', dclass: 'button_toggle', onpress : function () { $("#flex1").flexColToggler(this); } },  
        {separator: true},
        {name: 'Reset', dclass: 'button_reset', onpress : function () { $("#flex1").flexReset(); window.location.href = window.location.href; } },  
        {separator: true},
        {name: 'Edit', dclass: 'button_edit', onpress : editPress },
        {separator: true},
        {name: "Run Registry Data", dclass: 'button_table',  onpress : tableToggle },
        {separator: true},
        {name: "View in DBS", onpress : function () {
          goLink("https://cmsweb.cern.ch/dbs_discovery/getRunsFromRange?ajax=0&userMode=user&release=Any&group=Any&tier=Any&dbsInst=cms_dbs_prod_global&primType=Any&primD=Any&minRun={RUN_NUMBER}&maxRun={RUN_NUMBER}");
        }}
      ]
    };

    // apply settings
    $("#flex1").flexigrid(fp);

    if ($.cookie("flex_timer")) {
      $.cookie("flex_timer", null);
      $("div.button_refresh").each(timerToggle);
    }

    if ($.cookie("flex_table_summary")) {
      $.cookie("flex_table_summary", null);
      $("div.button_table").each(tableToggle);
    }

    if ($.cookie("flex_multiselect")) {
      $.cookie("flex_multiselect", null);
      $("div.button_select").each(selectToggle);
    }

    $("#help").tooltip({
      delay: 0,
      fade: 250,
      showURL: false,
      bodyHandler: function() {
        return $("#helphtml").html();
      }
    });

    $("#roles").tooltip({
      delay: 0,
      fade: 250,
      showURL: false,
      bodyHandler: function() {
        return "<%= WebUtils.getCertInfo(request).replaceAll(", ", "<br/>").replaceAll("\"", "") %>";
      }
    });

    // init
    $("#flex1").flexReload();

    $("#dumpdatamenu").menu({
      hoverOpenDelay: 200 
    });

    $("#batchupdatemenu").menu({
      hoverOpenDelay: 200 
    });

    $("#batch_updater_progressbar").progressBar({ barImage: 'media/img/progressbg_red.gif', boxImage: 'media/img/progressbar.gif', showText: true});

  });

function dumpData(intpl, tpl, mime) {

  var url = "runregisterdata";
  if ($.cookie("flex_table_summary")) url = "runinfodata";
  url += "?format=xml";
  if (intpl != '') url += "&intpl=" + intpl;
  if (tpl != '')   url += "&tpl=" + tpl;
  if (mime != '')  url += "&mime=" + mime;

  if ($.cookie("flex_multiselect") == "true") {
    if($("#flex1 .trSelected").length == 0) {
      alert("In Multi Select mode at least one run must be selected for dumping.");
      return;
    }
    var runs = "";
    var selruns = $("#flex1 .trSelected td[field=RUN_NUMBER] div");
    $.each(selruns, function(i,c){ 
      runs += "(" + $(c).text() + ")?"; 
    } );
    url += "&qtype=RUN_NUMBER";
    url += "&query=^" + runs + "$";
  } else {
    if (total >= 1000) {
      alert("Too many data: " + total + " records. \n Please narrow search and try again.");
      return;
    }
    if ($.cookie("flex_qtype")) url += "&qtype=" + $.cookie("flex_qtype");
    if ($.cookie("flex_query")) url += "&query=" + $.cookie("flex_query");
  }

  if ($.cookie("flex_sortname")) url += "&sortname=" + $.cookie("flex_sortname");
  if ($.cookie("flex_sortorder")) url += "&sortorder=" + $.cookie("flex_sortorder");
  window.open(url);
}

function changeStatusTo(status) {

  var len = $("#flex1 .trSelected").length;
  if(len == 0) {
    alert("At least one run must be selected for batch update.");
    return;
  }

  var tag = '';
  var status_name = status;
  if (status == 'COMPLETED') {
    tag = getTag();
    status_name += "(" + tag + ")";
  }

  if (!confirm("Do you wish to move to " + status_name + " selected " + $("#flex1 .trSelected").length + " run(s)?")) return;

  var cancel = false;
  var pbvalue = 0;
  $("#batchupdatemenu").hide();
  $(".menu-ul innerbox").hide();
  $("#batch_updater").show();
  $("#batch_updater_progressbar").progressBar(0);

  $.each($("#flex1 .trSelected td[field=RUN_NUMBER] div"), function(i, c){ 
    if (cancel) return;
    var run = $(c).text();
    $.ajax({
      type: "GET",
      async: false,
      url: "beditprovider?run_number=" + run + "&action=changestatus&status=" + status + "&tag=" + tag,
      processData: false,
      dataType: "xml",
      error: function(o) {
        if (!cancel) {
          if(!confirm("Error!\n\n" + o.responseText + "\n\nDo you wish to continue?")) {
            cancel = 1;
          }
        }
      },
      complete: function(ret) {
        pbvalue += 100 / len;
        $("#batch_updater_progressbar").progressBar(pbvalue);
      }
    });
  });

  $("#flex1").flexReload();
  $("#batch_updater").hide();
  $("#batchupdatemenu").show();
  $("#batch_updater_progressbar").progressBar(0);
}

</script>
</head>

<body>
  <table width="100%">
    <tr>
      <td id="logo_img"></td>
      <td id="logo">CMS DQM Run Registry</td>
      <td id="login">

<% if (WebUtils.hasLoggedRole(request, WebUtils.EXPERT)) { %>

        <span id="batch_updater" style="display: none;" >
          Updating: 
          <span id="batch_updater_progressbar"></span>
        </span>

        <span id="batchupdatemenu"><a href="#">Update Selected</a>
          <ul>
            <li><a href="#" onclick="changeStatusTo('ONLINE')">To ONLINE</a></li>
            <li><a href="#" onclick="changeStatusTo('OFFLINE')">To OFFLINE</a></li>
            <li><a href="#" onclick="changeStatusTo('SIGNOFF')">To SIGNOFF</a></li>
            <li><a href="#" onclick="changeStatusTo('COMPLETED')">To COMPLETED</a></li>
          </ul>
        </span>

        &nbsp;|&nbsp;

<% } %>        

        <span id="dumpdatamenu"><a href="#">Dump Data</a>
          <ul>
            <li><a href="#" onclick="dumpData('xml,elog', '', 'text/plain')">ELOG</a></li>
            <li><a href="#" onclick="dumpData('xml,twiki', '', 'text/plain')">TWiki</a></li>
            <li><a href="#" onclick="dumpData('xml,text_csv', '', 'text/plain')">Text&nbsp;(CSV)</a></li>
            <li><a href="#" onclick="dumpData('xml,text_tsv', '', 'text/plain')">Text&nbsp;(TSV)</a></li>
            <li><a href="#" onclick="dumpData('xml', '', 'text/xml')">XML</a></li>
            <li><a href="#" onclick="dumpData('xml,table', '', 'text/html')">HTML&nbsp;Table</a></li>
          </ul>
        </span>

        &nbsp;|&nbsp;

        <a id="help" href="help.html" target="_blank">Quick Help</a>

        &nbsp;|&nbsp;

	<a id="tutorial" 
href="https://twiki.cern.ch/twiki/bin/view/CMS/CMSDQMRunRegistry" 
target="_blank">Tutorial</a>

        &nbsp;|&nbsp;

        <% if (logged_in != null) { %>
          Logged in as <%= logged_in %> (<span id="roles"><%= logged_roles %></span>) - 
          <a href="#" onclick="logoutUser()">Logout</a>
        <% } else { %>
          <a href="#" onclick="location.href=location.href.replace(/^http:\/\//,'https:\/\/')">Login</a>
        <% } %>

      </td>
    </tr>
  </table>

  <jsp:include page="edit.jsp" />

  <table id="flex1"></table>

  <br/>
  <div>
    <b>Status reference: </b>
    <span class="statGOOD">&nbsp;&nbsp;&nbsp;&nbsp; - GOOD</span>
    <span class="statBAD">&nbsp;&nbsp;&nbsp;&nbsp; - BAD</span>
    <span class="statNOTSET">&nbsp;&nbsp;&nbsp;&nbsp; - NOTSET (offline subsystems only)</span>
    <span class="statEXCL">&nbsp;&nbsp;&nbsp;&nbsp; - EXCLUDED (online subsystems only)</span>
  </div>

  <iframe name="logout" width="1" height="1" src="" style="display:none;"></iframe>

  <div id="helphtml" style="display: none;">
    <jsp:include page="help.html" />
  </div>

<% if (logged_in != null) { %>
  <script language="javascript">
   setInterval("checkSession();", 20000);
  </script>
<% } %>

</body>
</html>

<!-- 
  NICE User:   <%= logged_in %>
  Certificate: <%= WebUtils.getCertInfo(request) %> 
  Roles:       <%= logged_roles %>
  Address:     <%= request.getRemoteAddr()%>
-->

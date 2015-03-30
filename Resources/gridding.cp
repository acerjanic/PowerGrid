<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='utf-8'>
<meta content='GitLab Community Edition' name='description'>
<title>
GPU_recon/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp at master - MRFIL  / GPU_recon | 
GitLab
</title>
<link href="/assets/favicon-baaa14bade1248aa6165e9d34e7d83c0.ico" rel="shortcut icon" type="image/vnd.microsoft.icon" />
<link href="/assets/application-2ac58be704fef592dcfa8a14124e2688.css" media="all" rel="stylesheet" />
<link href="/assets/print-1df3ea9b8ff148a6745321899e0cb213.css" media="print" rel="stylesheet" />
<script src="/assets/application-5f3c67da81199dd3da676a30746cf17f.js"></script>
<meta content="authenticity_token" name="csrf-param" />
<meta content="QMk16ncmfORs6tMHv3dq6t7rr22UsMHLdG6E0bli818=" name="csrf-token" />
<script type="text/javascript">
//<![CDATA[
window.gon={};gon.default_issues_tracker="gitlab";gon.api_version="v3";gon.relative_url_root="";gon.default_avatar_url="http://bioe-mrfil-07.bioen.illinois.edu/assets/no_avatar-adffbfe10d45b20495cd2a9b88974150.png";gon.current_user_id=3;gon.api_token="Ey9CCcQiJw5FBBQKsyGp";
//]]>
</script>
<meta content='width=device-width, initial-scale=1.0' name='viewport'>


</head>

<body class='ui_mars dark_theme project' data-page='projects:blob:show' data-project-id='46'>

<header class='navbar navbar-fixed-top navbar-gitlab'>
<div class='navbar-inner'>
<div class='container'>
<div class='app_logo'>
<a class="home has_bottom_tooltip" href="/" title="Dashboard"><h1>GITLAB</h1>
</a></div>
<h1 class='title'><span><a href="/groups/mrfil">MRFIL </a> / <a href="/mrfil/GPU_recon">GPU_recon</a></span></h1>
<button class='navbar-toggle' data-target='.navbar-collapse' data-toggle='collapse' type='button'>
<span class='sr-only'>Toggle navigation</span>
<i class='fa fa-bars'></i>
</button>
<div class='navbar-collapse collapse'>
<ul class='nav navbar-nav'>
<li class='hidden-sm hidden-xs'>
<div class='search'>
<form accept-charset="UTF-8" action="/search" class="navbar-form pull-left" method="get"><div style="display:none"><input name="utf8" type="hidden" value="&#x2713;" /></div>
<input class="search-input" id="search" name="search" placeholder="Search in this project" type="search" />
<input id="group_id" name="group_id" type="hidden" />
<input id="project_id" name="project_id" type="hidden" value="46" />
<input id="search_code" name="search_code" type="hidden" value="true" />
<input id="repository_ref" name="repository_ref" type="hidden" value="master" />

<div class='search-autocomplete-opts hide' data-autocomplete-path='/search/autocomplete' data-autocomplete-project-id='46' data-autocomplete-project-ref='master'></div>
</form>

</div>
<script>
  $('.search-input').on('keyup', function(e) {
    if (e.keyCode == 27) {
      $('.search-input').blur()
    }
  })
</script>

</li>
<li class='visible-sm visible-xs'>
<a class="has_bottom_tooltip" data-original-title="Search area" href="/search" title="Search"><i class='fa fa-search'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-original-title="Help" href="/help" title="Help"><i class='fa fa-question-circle'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-original-title="Public area" href="/explore" title="Explore"><i class='fa fa-globe'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-original-title="My snippets" href="/s/acerja2" title="My snippets"><i class='fa fa-clipboard'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-original-title="Admin area" href="/admin" title="Admin area"><i class='fa fa-cogs'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-original-title="New project" href="/projects/new" title="New project"><i class='fa fa-plus'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-original-title="Profile settings&quot;" href="/profile" title="Profile settings"><i class='fa fa-user'></i>
</a></li>
<li>
<a class="has_bottom_tooltip" data-method="delete" data-original-title="Logout" href="/users/sign_out" rel="nofollow" title="Logout"><i class='fa fa-sign-out'></i>
</a></li>
<li class='hidden-xs'>
<a class="profile-pic" href="/u/acerja2" id="profile-pic"><img alt="User activity" src="http://bioe-mrfil-07.bioen.illinois.edu//uploads/user/avatar/3/mrfil7.jpeg" />
</a></li>
</ul>
</div>
</div>
</div>
</header>


<script>
  GitLab.GfmAutoComplete.dataSource = "/mrfil/GPU_recon/autocomplete_sources?type=NilClass&type_id=master%2FIMPATIENT%2F3p1a%2FQuadro_FX_5600%2Ftoeplitz%2Fgridding.cpp"
  GitLab.GfmAutoComplete.setup();
</script>

<div class='page-with-sidebar'>
<div class='sidebar-wrapper'>
<ul class='project-navigation nav nav-sidebar'>
<li class="home"><a class="shortcuts-project" href="/mrfil/GPU_recon" title="Project"><i class='fa fa-dashboard'></i>
<span>
Project
</span>
</a></li><li class="active"><a class="shortcuts-tree" href="/mrfil/GPU_recon/tree/master"><i class='fa fa-files-o'></i>
<span>
Files
</span>
</a></li><li class=""><a class="shortcuts-commits" href="/mrfil/GPU_recon/commits/master"><i class='fa fa-history'></i>
<span>
Commits
</span>
</a></li><li class=""><a class="shortcuts-network" href="/mrfil/GPU_recon/network/master"><i class='fa fa-code-fork'></i>
<span>
Network
</span>
</a></li><li class=""><a class="shortcuts-graphs" href="/mrfil/GPU_recon/graphs/master"><i class='fa fa-area-chart'></i>
<span>
Graphs
</span>
</a></li><li class=""><a class="shortcuts-issues" href="/mrfil/GPU_recon/issues"><i class='fa fa-exclamation-circle'></i>
<span>
Issues
<span class='count issue_counter'>0</span>
</span>
</a></li><li class=""><a class="shortcuts-merge_requests" href="/mrfil/GPU_recon/merge_requests"><i class='fa fa-tasks'></i>
<span>
Merge Requests
<span class='count merge_counter'>0</span>
</span>
</a></li><li class=""><a class="shortcuts-wiki" href="/mrfil/GPU_recon/wikis/home"><i class='fa fa-book'></i>
<span>
Wiki
</span>
</a></li><li class="separate-item"><a class="stat-tab tab no-highlight" href="/mrfil/GPU_recon/edit"><i class='fa fa-cogs'></i>
<span>
Settings
<i class='fa fa-angle-down'></i>
</span>
</a></li></ul>

</div>
<div class='content-wrapper'>
<div class='container-fluid'>
<div class='content'>
<div class='flash-container'>
</div>

<div class='clearfix'>
<div class='tree-ref-holder'>
<form accept-charset="UTF-8" action="/mrfil/GPU_recon/refs/switch" class="project-refs-form" method="get"><div style="display:none"><input name="utf8" type="hidden" value="&#x2713;" /></div>
<select class="project-refs-select select2 select2-sm" id="ref" name="ref"><optgroup label="Branches"><option selected="selected" value="master">master</option></optgroup><optgroup label="Tags"></optgroup></select>
<input id="destination" name="destination" type="hidden" value="blob" />
<input id="path" name="path" type="hidden" value="IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp" />
</form>


</div>
<div class='tree-holder' id='tree-holder'>
<ul class='breadcrumb repo-breadcrumb'>
<li>
<i class='fa fa-angle-right'></i>
<a href="/mrfil/GPU_recon/tree/master">GPU_recon
</a></li>
<li>
<a href="/mrfil/GPU_recon/tree/master/IMPATIENT">IMPATIENT</a>
</li>
<li>
<a href="/mrfil/GPU_recon/tree/master/IMPATIENT/3p1a">3p1a</a>
</li>
<li>
<a href="/mrfil/GPU_recon/tree/master/IMPATIENT/3p1a/Quadro_FX_5600">Quadro_FX_5600</a>
</li>
<li>
<a href="/mrfil/GPU_recon/tree/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz">toeplitz</a>
</li>
<li>
<a href="/mrfil/GPU_recon/blob/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp"><strong>
gridding.cpp
</strong>
</a></li>
</ul>
<ul class='blob-commit-info bs-callout bs-callout-info hidden-xs'>
<li class='commit js-toggle-container'>
<div class='commit-row-title'>
<a class="commit_short_id" href="/mrfil/GPU_recon/commit/e0f8cdb8273f1aeda60c70364069d47c14eedae8">e0f8cdb8</a>
&nbsp;
<span class='str-truncated'>
<a class="commit-row-message" href="/mrfil/GPU_recon/commit/e0f8cdb8273f1aeda60c70364069d47c14eedae8">Initial Commit</a>
</span>
<a class="pull-right" href="/mrfil/GPU_recon/tree/e0f8cdb8273f1aeda60c70364069d47c14eedae8">Browse Code »</a>
<div class='notes_count'>
</div>
</div>
<div class='commit-row-info'>
<a class="commit-author-link has_tooltip" data-original-title="holtrop1@illinois.edu" href="/u/holtrop1"><img alt="" class="avatar s16" src="http://bioe-mrfil-07.bioen.illinois.edu//uploads/user/avatar/2/mrfil9.jpeg" width="16" /> <span class="commit-author-name">Joe Holtrop</span></a>
<div class='committed_ago'>
<time class='time_ago' data-placement='top' data-toggle='tooltip' datetime='2015-03-19T23:24:22Z' title='Mar 19, 2015 6:24pm'>2015-03-19 18:24:22 -0500</time>
<script>$('.time_ago').timeago().tooltip()</script>
 &nbsp;
</div>
</div>
</li>

</ul>
<div class='tree-content-holder' id='tree-content-holder'>
<article class='file-holder'>
<div class='file-title clearfix'>
<i class='fa fa-file'></i>
<span class='file_name'>
gridding.cpp
<small>5.36 KB</small>
</span>
<span class='options hidden-xs'><div class='btn-group tree-btn-group'>
<a class="btn btn-small" href="/mrfil/GPU_recon/edit/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp">Edit</a>
<a class="btn btn-small" href="/mrfil/GPU_recon/raw/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp" target="_blank">Raw</a>
<a class="btn btn-small" href="/mrfil/GPU_recon/blame/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp">Blame</a>
<a class="btn btn-small" href="/mrfil/GPU_recon/commits/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp">History</a>
<a class="btn btn-small" href="/mrfil/GPU_recon/blob/e0f8cdb8273f1aeda60c70364069d47c14eedae8/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp">Permalink</a>
</div>
<button class="remove-blob btn btn-small btn-remove" data-target="#modal-remove-blob" data-toggle="modal" name="button" type="submit">Remove
</button></span>
</div>
<div class='file-content code'>
<div class='dark highlighted-data'>
<div class='line-numbers'>
<a href="#L1" id="L1" rel="#L1"><i class='fa fa-link'></i>
1
</a><a href="#L2" id="L2" rel="#L2"><i class='fa fa-link'></i>
2
</a><a href="#L3" id="L3" rel="#L3"><i class='fa fa-link'></i>
3
</a><a href="#L4" id="L4" rel="#L4"><i class='fa fa-link'></i>
4
</a><a href="#L5" id="L5" rel="#L5"><i class='fa fa-link'></i>
5
</a><a href="#L6" id="L6" rel="#L6"><i class='fa fa-link'></i>
6
</a><a href="#L7" id="L7" rel="#L7"><i class='fa fa-link'></i>
7
</a><a href="#L8" id="L8" rel="#L8"><i class='fa fa-link'></i>
8
</a><a href="#L9" id="L9" rel="#L9"><i class='fa fa-link'></i>
9
</a><a href="#L10" id="L10" rel="#L10"><i class='fa fa-link'></i>
10
</a><a href="#L11" id="L11" rel="#L11"><i class='fa fa-link'></i>
11
</a><a href="#L12" id="L12" rel="#L12"><i class='fa fa-link'></i>
12
</a><a href="#L13" id="L13" rel="#L13"><i class='fa fa-link'></i>
13
</a><a href="#L14" id="L14" rel="#L14"><i class='fa fa-link'></i>
14
</a><a href="#L15" id="L15" rel="#L15"><i class='fa fa-link'></i>
15
</a><a href="#L16" id="L16" rel="#L16"><i class='fa fa-link'></i>
16
</a><a href="#L17" id="L17" rel="#L17"><i class='fa fa-link'></i>
17
</a><a href="#L18" id="L18" rel="#L18"><i class='fa fa-link'></i>
18
</a><a href="#L19" id="L19" rel="#L19"><i class='fa fa-link'></i>
19
</a><a href="#L20" id="L20" rel="#L20"><i class='fa fa-link'></i>
20
</a><a href="#L21" id="L21" rel="#L21"><i class='fa fa-link'></i>
21
</a><a href="#L22" id="L22" rel="#L22"><i class='fa fa-link'></i>
22
</a><a href="#L23" id="L23" rel="#L23"><i class='fa fa-link'></i>
23
</a><a href="#L24" id="L24" rel="#L24"><i class='fa fa-link'></i>
24
</a><a href="#L25" id="L25" rel="#L25"><i class='fa fa-link'></i>
25
</a><a href="#L26" id="L26" rel="#L26"><i class='fa fa-link'></i>
26
</a><a href="#L27" id="L27" rel="#L27"><i class='fa fa-link'></i>
27
</a><a href="#L28" id="L28" rel="#L28"><i class='fa fa-link'></i>
28
</a><a href="#L29" id="L29" rel="#L29"><i class='fa fa-link'></i>
29
</a><a href="#L30" id="L30" rel="#L30"><i class='fa fa-link'></i>
30
</a><a href="#L31" id="L31" rel="#L31"><i class='fa fa-link'></i>
31
</a><a href="#L32" id="L32" rel="#L32"><i class='fa fa-link'></i>
32
</a><a href="#L33" id="L33" rel="#L33"><i class='fa fa-link'></i>
33
</a><a href="#L34" id="L34" rel="#L34"><i class='fa fa-link'></i>
34
</a><a href="#L35" id="L35" rel="#L35"><i class='fa fa-link'></i>
35
</a><a href="#L36" id="L36" rel="#L36"><i class='fa fa-link'></i>
36
</a><a href="#L37" id="L37" rel="#L37"><i class='fa fa-link'></i>
37
</a><a href="#L38" id="L38" rel="#L38"><i class='fa fa-link'></i>
38
</a><a href="#L39" id="L39" rel="#L39"><i class='fa fa-link'></i>
39
</a><a href="#L40" id="L40" rel="#L40"><i class='fa fa-link'></i>
40
</a><a href="#L41" id="L41" rel="#L41"><i class='fa fa-link'></i>
41
</a><a href="#L42" id="L42" rel="#L42"><i class='fa fa-link'></i>
42
</a><a href="#L43" id="L43" rel="#L43"><i class='fa fa-link'></i>
43
</a><a href="#L44" id="L44" rel="#L44"><i class='fa fa-link'></i>
44
</a><a href="#L45" id="L45" rel="#L45"><i class='fa fa-link'></i>
45
</a><a href="#L46" id="L46" rel="#L46"><i class='fa fa-link'></i>
46
</a><a href="#L47" id="L47" rel="#L47"><i class='fa fa-link'></i>
47
</a><a href="#L48" id="L48" rel="#L48"><i class='fa fa-link'></i>
48
</a><a href="#L49" id="L49" rel="#L49"><i class='fa fa-link'></i>
49
</a><a href="#L50" id="L50" rel="#L50"><i class='fa fa-link'></i>
50
</a><a href="#L51" id="L51" rel="#L51"><i class='fa fa-link'></i>
51
</a><a href="#L52" id="L52" rel="#L52"><i class='fa fa-link'></i>
52
</a><a href="#L53" id="L53" rel="#L53"><i class='fa fa-link'></i>
53
</a><a href="#L54" id="L54" rel="#L54"><i class='fa fa-link'></i>
54
</a><a href="#L55" id="L55" rel="#L55"><i class='fa fa-link'></i>
55
</a><a href="#L56" id="L56" rel="#L56"><i class='fa fa-link'></i>
56
</a><a href="#L57" id="L57" rel="#L57"><i class='fa fa-link'></i>
57
</a><a href="#L58" id="L58" rel="#L58"><i class='fa fa-link'></i>
58
</a><a href="#L59" id="L59" rel="#L59"><i class='fa fa-link'></i>
59
</a><a href="#L60" id="L60" rel="#L60"><i class='fa fa-link'></i>
60
</a><a href="#L61" id="L61" rel="#L61"><i class='fa fa-link'></i>
61
</a><a href="#L62" id="L62" rel="#L62"><i class='fa fa-link'></i>
62
</a><a href="#L63" id="L63" rel="#L63"><i class='fa fa-link'></i>
63
</a><a href="#L64" id="L64" rel="#L64"><i class='fa fa-link'></i>
64
</a><a href="#L65" id="L65" rel="#L65"><i class='fa fa-link'></i>
65
</a><a href="#L66" id="L66" rel="#L66"><i class='fa fa-link'></i>
66
</a><a href="#L67" id="L67" rel="#L67"><i class='fa fa-link'></i>
67
</a><a href="#L68" id="L68" rel="#L68"><i class='fa fa-link'></i>
68
</a><a href="#L69" id="L69" rel="#L69"><i class='fa fa-link'></i>
69
</a><a href="#L70" id="L70" rel="#L70"><i class='fa fa-link'></i>
70
</a><a href="#L71" id="L71" rel="#L71"><i class='fa fa-link'></i>
71
</a><a href="#L72" id="L72" rel="#L72"><i class='fa fa-link'></i>
72
</a><a href="#L73" id="L73" rel="#L73"><i class='fa fa-link'></i>
73
</a><a href="#L74" id="L74" rel="#L74"><i class='fa fa-link'></i>
74
</a><a href="#L75" id="L75" rel="#L75"><i class='fa fa-link'></i>
75
</a><a href="#L76" id="L76" rel="#L76"><i class='fa fa-link'></i>
76
</a><a href="#L77" id="L77" rel="#L77"><i class='fa fa-link'></i>
77
</a><a href="#L78" id="L78" rel="#L78"><i class='fa fa-link'></i>
78
</a><a href="#L79" id="L79" rel="#L79"><i class='fa fa-link'></i>
79
</a><a href="#L80" id="L80" rel="#L80"><i class='fa fa-link'></i>
80
</a><a href="#L81" id="L81" rel="#L81"><i class='fa fa-link'></i>
81
</a><a href="#L82" id="L82" rel="#L82"><i class='fa fa-link'></i>
82
</a><a href="#L83" id="L83" rel="#L83"><i class='fa fa-link'></i>
83
</a><a href="#L84" id="L84" rel="#L84"><i class='fa fa-link'></i>
84
</a><a href="#L85" id="L85" rel="#L85"><i class='fa fa-link'></i>
85
</a><a href="#L86" id="L86" rel="#L86"><i class='fa fa-link'></i>
86
</a><a href="#L87" id="L87" rel="#L87"><i class='fa fa-link'></i>
87
</a><a href="#L88" id="L88" rel="#L88"><i class='fa fa-link'></i>
88
</a><a href="#L89" id="L89" rel="#L89"><i class='fa fa-link'></i>
89
</a><a href="#L90" id="L90" rel="#L90"><i class='fa fa-link'></i>
90
</a><a href="#L91" id="L91" rel="#L91"><i class='fa fa-link'></i>
91
</a><a href="#L92" id="L92" rel="#L92"><i class='fa fa-link'></i>
92
</a><a href="#L93" id="L93" rel="#L93"><i class='fa fa-link'></i>
93
</a><a href="#L94" id="L94" rel="#L94"><i class='fa fa-link'></i>
94
</a><a href="#L95" id="L95" rel="#L95"><i class='fa fa-link'></i>
95
</a><a href="#L96" id="L96" rel="#L96"><i class='fa fa-link'></i>
96
</a><a href="#L97" id="L97" rel="#L97"><i class='fa fa-link'></i>
97
</a><a href="#L98" id="L98" rel="#L98"><i class='fa fa-link'></i>
98
</a><a href="#L99" id="L99" rel="#L99"><i class='fa fa-link'></i>
99
</a><a href="#L100" id="L100" rel="#L100"><i class='fa fa-link'></i>
100
</a><a href="#L101" id="L101" rel="#L101"><i class='fa fa-link'></i>
101
</a><a href="#L102" id="L102" rel="#L102"><i class='fa fa-link'></i>
102
</a><a href="#L103" id="L103" rel="#L103"><i class='fa fa-link'></i>
103
</a><a href="#L104" id="L104" rel="#L104"><i class='fa fa-link'></i>
104
</a><a href="#L105" id="L105" rel="#L105"><i class='fa fa-link'></i>
105
</a><a href="#L106" id="L106" rel="#L106"><i class='fa fa-link'></i>
106
</a><a href="#L107" id="L107" rel="#L107"><i class='fa fa-link'></i>
107
</a><a href="#L108" id="L108" rel="#L108"><i class='fa fa-link'></i>
108
</a><a href="#L109" id="L109" rel="#L109"><i class='fa fa-link'></i>
109
</a><a href="#L110" id="L110" rel="#L110"><i class='fa fa-link'></i>
110
</a><a href="#L111" id="L111" rel="#L111"><i class='fa fa-link'></i>
111
</a><a href="#L112" id="L112" rel="#L112"><i class='fa fa-link'></i>
112
</a><a href="#L113" id="L113" rel="#L113"><i class='fa fa-link'></i>
113
</a><a href="#L114" id="L114" rel="#L114"><i class='fa fa-link'></i>
114
</a><a href="#L115" id="L115" rel="#L115"><i class='fa fa-link'></i>
115
</a><a href="#L116" id="L116" rel="#L116"><i class='fa fa-link'></i>
116
</a><a href="#L117" id="L117" rel="#L117"><i class='fa fa-link'></i>
117
</a><a href="#L118" id="L118" rel="#L118"><i class='fa fa-link'></i>
118
</a><a href="#L119" id="L119" rel="#L119"><i class='fa fa-link'></i>
119
</a><a href="#L120" id="L120" rel="#L120"><i class='fa fa-link'></i>
120
</a><a href="#L121" id="L121" rel="#L121"><i class='fa fa-link'></i>
121
</a><a href="#L122" id="L122" rel="#L122"><i class='fa fa-link'></i>
122
</a><a href="#L123" id="L123" rel="#L123"><i class='fa fa-link'></i>
123
</a><a href="#L124" id="L124" rel="#L124"><i class='fa fa-link'></i>
124
</a><a href="#L125" id="L125" rel="#L125"><i class='fa fa-link'></i>
125
</a><a href="#L126" id="L126" rel="#L126"><i class='fa fa-link'></i>
126
</a><a href="#L127" id="L127" rel="#L127"><i class='fa fa-link'></i>
127
</a><a href="#L128" id="L128" rel="#L128"><i class='fa fa-link'></i>
128
</a><a href="#L129" id="L129" rel="#L129"><i class='fa fa-link'></i>
129
</a><a href="#L130" id="L130" rel="#L130"><i class='fa fa-link'></i>
130
</a><a href="#L131" id="L131" rel="#L131"><i class='fa fa-link'></i>
131
</a><a href="#L132" id="L132" rel="#L132"><i class='fa fa-link'></i>
132
</a><a href="#L133" id="L133" rel="#L133"><i class='fa fa-link'></i>
133
</a><a href="#L134" id="L134" rel="#L134"><i class='fa fa-link'></i>
134
</a><a href="#L135" id="L135" rel="#L135"><i class='fa fa-link'></i>
135
</a><a href="#L136" id="L136" rel="#L136"><i class='fa fa-link'></i>
136
</a><a href="#L137" id="L137" rel="#L137"><i class='fa fa-link'></i>
137
</a><a href="#L138" id="L138" rel="#L138"><i class='fa fa-link'></i>
138
</a><a href="#L139" id="L139" rel="#L139"><i class='fa fa-link'></i>
139
</a><a href="#L140" id="L140" rel="#L140"><i class='fa fa-link'></i>
140
</a><a href="#L141" id="L141" rel="#L141"><i class='fa fa-link'></i>
141
</a><a href="#L142" id="L142" rel="#L142"><i class='fa fa-link'></i>
142
</a><a href="#L143" id="L143" rel="#L143"><i class='fa fa-link'></i>
143
</a><a href="#L144" id="L144" rel="#L144"><i class='fa fa-link'></i>
144
</a><a href="#L145" id="L145" rel="#L145"><i class='fa fa-link'></i>
145
</a><a href="#L146" id="L146" rel="#L146"><i class='fa fa-link'></i>
146
</a><a href="#L147" id="L147" rel="#L147"><i class='fa fa-link'></i>
147
</a><a href="#L148" id="L148" rel="#L148"><i class='fa fa-link'></i>
148
</a><a href="#L149" id="L149" rel="#L149"><i class='fa fa-link'></i>
149
</a><a href="#L150" id="L150" rel="#L150"><i class='fa fa-link'></i>
150
</a><a href="#L151" id="L151" rel="#L151"><i class='fa fa-link'></i>
151
</a><a href="#L152" id="L152" rel="#L152"><i class='fa fa-link'></i>
152
</a><a href="#L153" id="L153" rel="#L153"><i class='fa fa-link'></i>
153
</a><a href="#L154" id="L154" rel="#L154"><i class='fa fa-link'></i>
154
</a><a href="#L155" id="L155" rel="#L155"><i class='fa fa-link'></i>
155
</a><a href="#L156" id="L156" rel="#L156"><i class='fa fa-link'></i>
156
</a><a href="#L157" id="L157" rel="#L157"><i class='fa fa-link'></i>
157
</a><a href="#L158" id="L158" rel="#L158"><i class='fa fa-link'></i>
158
</a><a href="#L159" id="L159" rel="#L159"><i class='fa fa-link'></i>
159
</a><a href="#L160" id="L160" rel="#L160"><i class='fa fa-link'></i>
160
</a><a href="#L161" id="L161" rel="#L161"><i class='fa fa-link'></i>
161
</a><a href="#L162" id="L162" rel="#L162"><i class='fa fa-link'></i>
162
</a><a href="#L163" id="L163" rel="#L163"><i class='fa fa-link'></i>
163
</a><a href="#L164" id="L164" rel="#L164"><i class='fa fa-link'></i>
164
</a><a href="#L165" id="L165" rel="#L165"><i class='fa fa-link'></i>
165
</a><a href="#L166" id="L166" rel="#L166"><i class='fa fa-link'></i>
166
</a><a href="#L167" id="L167" rel="#L167"><i class='fa fa-link'></i>
167
</a><a href="#L168" id="L168" rel="#L168"><i class='fa fa-link'></i>
168
</a><a href="#L169" id="L169" rel="#L169"><i class='fa fa-link'></i>
169
</a><a href="#L170" id="L170" rel="#L170"><i class='fa fa-link'></i>
170
</a><a href="#L171" id="L171" rel="#L171"><i class='fa fa-link'></i>
171
</a><a href="#L172" id="L172" rel="#L172"><i class='fa fa-link'></i>
172
</a><a href="#L173" id="L173" rel="#L173"><i class='fa fa-link'></i>
173
</a><a href="#L174" id="L174" rel="#L174"><i class='fa fa-link'></i>
174
</a><a href="#L175" id="L175" rel="#L175"><i class='fa fa-link'></i>
175
</a><a href="#L176" id="L176" rel="#L176"><i class='fa fa-link'></i>
176
</a><a href="#L177" id="L177" rel="#L177"><i class='fa fa-link'></i>
177
</a><a href="#L178" id="L178" rel="#L178"><i class='fa fa-link'></i>
178
</a><a href="#L179" id="L179" rel="#L179"><i class='fa fa-link'></i>
179
</a><a href="#L180" id="L180" rel="#L180"><i class='fa fa-link'></i>
180
</a><a href="#L181" id="L181" rel="#L181"><i class='fa fa-link'></i>
181
</a><a href="#L182" id="L182" rel="#L182"><i class='fa fa-link'></i>
182
</a><a href="#L183" id="L183" rel="#L183"><i class='fa fa-link'></i>
183
</a><a href="#L184" id="L184" rel="#L184"><i class='fa fa-link'></i>
184
</a><a href="#L185" id="L185" rel="#L185"><i class='fa fa-link'></i>
185
</a><a href="#L186" id="L186" rel="#L186"><i class='fa fa-link'></i>
186
</a><a href="#L187" id="L187" rel="#L187"><i class='fa fa-link'></i>
187
</a><a href="#L188" id="L188" rel="#L188"><i class='fa fa-link'></i>
188
</a><a href="#L189" id="L189" rel="#L189"><i class='fa fa-link'></i>
189
</a><a href="#L190" id="L190" rel="#L190"><i class='fa fa-link'></i>
190
</a><a href="#L191" id="L191" rel="#L191"><i class='fa fa-link'></i>
191
</a><a href="#L192" id="L192" rel="#L192"><i class='fa fa-link'></i>
192
</a><a href="#L193" id="L193" rel="#L193"><i class='fa fa-link'></i>
193
</a><a href="#L194" id="L194" rel="#L194"><i class='fa fa-link'></i>
194
</a><a href="#L195" id="L195" rel="#L195"><i class='fa fa-link'></i>
195
</a><a href="#L196" id="L196" rel="#L196"><i class='fa fa-link'></i>
196
</a><a href="#L197" id="L197" rel="#L197"><i class='fa fa-link'></i>
197
</a><a href="#L198" id="L198" rel="#L198"><i class='fa fa-link'></i>
198
</a><a href="#L199" id="L199" rel="#L199"><i class='fa fa-link'></i>
199
</a><a href="#L200" id="L200" rel="#L200"><i class='fa fa-link'></i>
200
</a><a href="#L201" id="L201" rel="#L201"><i class='fa fa-link'></i>
201
</a><a href="#L202" id="L202" rel="#L202"><i class='fa fa-link'></i>
202
</a><a href="#L203" id="L203" rel="#L203"><i class='fa fa-link'></i>
203
</a><a href="#L204" id="L204" rel="#L204"><i class='fa fa-link'></i>
204
</a><a href="#L205" id="L205" rel="#L205"><i class='fa fa-link'></i>
205
</a><a href="#L206" id="L206" rel="#L206"><i class='fa fa-link'></i>
206
</a><a href="#L207" id="L207" rel="#L207"><i class='fa fa-link'></i>
207
</a><a href="#L208" id="L208" rel="#L208"><i class='fa fa-link'></i>
208
</a><a href="#L209" id="L209" rel="#L209"><i class='fa fa-link'></i>
209
</a><a href="#L210" id="L210" rel="#L210"><i class='fa fa-link'></i>
210
</a><a href="#L211" id="L211" rel="#L211"><i class='fa fa-link'></i>
211
</a><a href="#L212" id="L212" rel="#L212"><i class='fa fa-link'></i>
212
</a></div>
<div class='highlight'>
<pre><code class='gridding.cpp'>#include &lt;gridding.h&gt;
#include &lt;math.h&gt;
#include &lt;utils.h&gt;
#include &lt;stdio.h&gt;
#include &lt;string.h&gt;
#include &lt;stdlib.h&gt;

// t0 is there b/c t.dat does not start with 0.0f.
static    __host__ __device__ 
float hanning_d(float tm, float tau, float l, float t0)
{
    float taul = tau * l;
    float result;
    if ( fabs(tm - taul - t0) &lt; tau ) {
        result = 0.5f + 0.5f * cosf(PI * (tm - taul - t0) / tau);
    } else {
        result = 0.0f;
    }
    //FIXME:
    //result = 1.0f;
    return result;
}

///* Jiading GAI - GRIDDING - BEGIN
/*From Numerical Recipes in C, 2nd Edition*/
static float bessi0(float x)
{
    float ax,ans;
    float y;
    
    if ((ax=fabs(x)) &lt; 3.75)
    {
        y=x/3.75;
        y=y*y;
        ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+
            y*(0.360768e-1+y*0.45813e-2)))));
    }
    else
    {
        y=3.75/ax;
        ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+
             y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+
             y*(-0.1647633e-1+y*0.392377e-2))))))));
    }
    return ans;
}

	void 
binning_kernel1_CPU(unsigned int n, ReconstructionSample *Sample_h,
                    unsigned int *numPts_h, parameters params)
{
  unsigned int gridNumElems = params.gridSize[0] * 
                              params.gridSize[1] * 
                              params.gridSize[2] ;
  unsigned int binsize = params.binsize;
  float gridOS = params.gridOS;
  unsigned int Nxy_grid = params.gridSize[0] * 
                          params.gridSize[1] ;

  for(unsigned int i=0;i&lt;n;i++)
  {

     ReconstructionSample pt;
     unsigned int binIdx;

     pt = Sample_h[i];
     pt.kX = (gridOS)*(pt.kX+((float)params.imageSize[0])/2.0f);
     pt.kY = (gridOS)*(pt.kY+((float)params.imageSize[1])/2.0f);
	 if(1==params.imageSize[2])
     {
	   pt.kZ = 0.0f;
	 }
	 else
     {
       pt.kZ = (gridOS)*(pt.kZ+((float)params.imageSize[2])/2.0f);
	 }

	 /*
	 // Clamp k trajectories between [0,gridSize]
	 // B/c not all data given to me are normalized 
	 // between [-imageSize/2,imageSize/2].
	 if( pt.kZ &lt; 0.0f )
		 pt.kZ = 0.0f;
	 if( pt.kZ &gt; (float)params.gridSize[2] )
		 pt.kZ = (float)params.gridSize[2];

	 if( pt.kX &lt; 0.0f )
		 pt.kX = 0.0f;
	 if( pt.kX &gt; (float)params.gridSize[0] )
		 pt.kX = (float)params.gridSize[0];

	 if( pt.kY &lt; 0.0f )
		 pt.kY = 0.0f;
	 if( pt.kY &gt; (float)params.gridSize[1] )
		 pt.kY = (float)params.gridSize[1];
     // */

     binIdx = (unsigned int)(pt.kZ)*Nxy_grid + 
              (unsigned int)(pt.kY)*params.gridSize[0] + 
              (unsigned int)(pt.kX);

     if(numPts_h[binIdx]&lt;binsize)
     {
        numPts_h[binIdx] += 1;
     }
  }  
}

  void
scanLargeArray_CPU(int len, unsigned int *input)
{
  unsigned int *output = (unsigned int*) malloc(len*sizeof(unsigned int));
  output[0] = 0;

  for(int i=1;i&lt;len;i++)
    output[i] = output[i-1] + input[i-1];

  memcpy(input, output, len*sizeof(unsigned int));
  free(output);
}

   void
binning_kernel2_CPU(unsigned int n, ReconstructionSample *Sample, unsigned int *binStartAddr,
                    unsigned int *numPts, parameters params, samplePtArray SampleArray,
                    int&amp; CPUbinSize, int *CPUbin)
{
   unsigned int gridNumElems = params.gridSize[0] * 
                               params.gridSize[1] * 
                               params.gridSize[2] ;
   unsigned int binsize = params.binsize;
   float gridOS = params.gridOS;
   unsigned int Nxy_grid = params.gridSize[0] * 
                           params.gridSize[1] ;

   for(unsigned int i=0;i&lt;n;i++)
   {
      ReconstructionSample pt;
      unsigned int binIdx;
      unsigned int cap;

      pt = Sample[i];
      pt.kX = (gridOS)*(pt.kX+((float)params.imageSize[0])/2.0f);
      pt.kY = (gridOS)*(pt.kY+((float)params.imageSize[1])/2.0f);
	  if(1==params.imageSize[2])
	  {
        pt.kZ = 0.0f;
	  }
	  else
	  {
        pt.kZ = (gridOS)*(pt.kZ+((float)params.imageSize[2])/2.0f);
	  }

	  /*
	  // Clamp k trajectories between [0,gridSize].
	  // B/c not all input data given to me are normalized 
	  // between [-imageSize/2,imageSize/2].
      if( pt.kZ &lt; 0.0f )
	      pt.kZ = 0.0f;
      if( pt.kZ &gt; (float)params.gridSize[2] )
	      pt.kZ = (float)params.gridSize[2];

      if( pt.kX &lt; 0.0f )
	      pt.kX = 0.0f;
      if( pt.kX &gt; (float)params.gridSize[0] )
	      pt.kX = (float)params.gridSize[0];

	 if( pt.kY &lt; 0.0f )
		 pt.kY = 0.0f;
	 if( pt.kY &gt; (float)params.gridSize[1] )
		 pt.kY = (float)params.gridSize[1];
	 // */

      binIdx = (unsigned int)(pt.kZ)*Nxy_grid + 
               (unsigned int)(pt.kY)*params.gridSize[0] + 
               (unsigned int)(pt.kX);

      if(numPts[binIdx]&lt;binsize)
      {
         cap = numPts[binIdx];
         numPts[binIdx] += 1;
         
         float2 data;
         data.x = pt.real;
         data.y = pt.imag;
 
         float2 loc1;
         loc1.x = pt.kX;
         loc1.y = pt.kY;
 
         float2 loc2;
         loc2.x = pt.kZ;
         loc2.y = pt.sdc;

         float2 loc3;
		 loc3.x = pt.t;
		 //loc3.y = pt.dummy;

         SampleArray.data[binStartAddr[binIdx]+cap] = data;
         SampleArray.loc1[binStartAddr[binIdx]+cap] = loc1;
         SampleArray.loc2[binStartAddr[binIdx]+cap] = loc2;
         SampleArray.loc3[binStartAddr[binIdx]+cap] = loc3;
      }
      else
      {
         cap = CPUbinSize;
         CPUbinSize += 1;
         CPUbin[cap] = i;
      }
   }  
}

// Jiading GAI - GRIDDING - END */
</code></pre>
</div>
</div>

</div>

</article>
</div>

</div>
<div class='modal hide' id='modal-remove-blob'>
<div class='modal-dialog'>
<div class='modal-content'>
<div class='modal-header'>
<a class='close' data-dismiss='modal' href='#'>×</a>
<h3 class='page-title'>Remove gridding.cpp</h3>
<p class='light'>
From branch
<strong>master</strong>
</p>
</div>
<div class='modal-body'>
<form accept-charset="UTF-8" action="/mrfil/GPU_recon/blob/master/IMPATIENT/3p1a/Quadro_FX_5600/toeplitz/gridding.cpp" class="form-horizontal" method="post"><div style="display:none"><input name="utf8" type="hidden" value="&#x2713;" /><input name="_method" type="hidden" value="delete" /><input name="authenticity_token" type="hidden" value="QMk16ncmfORs6tMHv3dq6t7rr22UsMHLdG6E0bli818=" /></div>
<div class='form-group commit_message-group'>
<label class="control-label" for="commit_message">Commit message
</label><div class='col-sm-10'>
<div class='commit-message-container'>
<div class='max-width-marker'></div>
<textarea class="form-control" id="commit_message" name="commit_message" placeholder="Removed this file because..." required="required" rows="3">
</textarea>
</div>
</div>
</div>

<div class='form-group'>
<div class='col-sm-2'></div>
<div class='col-sm-10'>
<button class="btn btn-remove btn-remove-file" name="button" type="submit">Remove file</button>
<a class="btn btn-cancel" data-dismiss="modal" href="#">Cancel</a>
</div>
</div>
</form>

</div>
</div>
</div>
</div>
<script>
  disableButtonIfEmptyField('#commit_message', '.btn-remove-file')
</script>


</div>
</div>
</div>
</div>
</div>


</body>
</html>

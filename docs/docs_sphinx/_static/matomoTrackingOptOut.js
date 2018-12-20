/******************************************************
 * Matomo Opt-out
 * 2018 Clive Beckett, https://coding.musikinsnetz.de/
 * licensed under the GNU GENERAL PUBLIC LICENSE, v3
 * Details and README under
 * https://github.com/clivebeckett/matomo-opt-out
 * Check README for the additional code you’ll need!
 ******************************************************/

/**
 * set whether you have your Matomo installation respect the browser’s DoNotTrack option
 * @var doNotTrackRespected boolean
 */
var doNotTrackRespected = true;

/**
 * language snippets
 * just duplicate one of the language blocks to add other languages
 * @object langSnippets
 */
var langSnippets = {
	"en": {
		"trackingActive": "Currently your visit to this website <strong>is being anonymously tracked</strong> by the Matomo Web Analytics Tool. Uncheck this box to stop being tracked.",
		"trackingInActive": "Currently your visit to this website <strong>is not being tracked</strong> by the Matomo Web Analytics Tool. Check this box to activate anonymous tracking and help us improve our website.",
		"doNotTrackActive": "You have activated the <em>Do Not Track</em> option in your browser settings. This setting is being respected by the Matomo Web Analytics Tool on our website. There’s no need for you to opt-out of this website’s data tracking.",
		"localStorageNotAvailable": "Your browser appears to be too old to support this feature. For more privacy control please update your browser."
	},
	"de": {
		"trackingActive": "Ihr Besuch auf dieser Website wird derzeit durch das Matomo-Webanalyse-Tool <strong>anonym erfasst</strong>. Klicken Sie hier um die Datenerfassung zu beenden.",
		"trackingInActive": "Ihr Besuch auf dieser Website wird derzeit <strong>nicht</strong> durch das Matomo-Webanalyse-Tool <strong>erfasst</strong>. Klicken Sie hier um anonyme Datenerfassung zu aktivieren und uns dabei zu helfen, unsere Website zu verbessern.",
		"doNotTrackActive": "Sie haben die <em>Do-Not-Track</em>-Option in Ihren Browser-Einstellungen aktiviert. Diese Einstellung wird von unserem Matomo-Webanalyse-Tool respektiert. Es ist daher nicht notwendig, die Datenerhebung für diese Website zu deaktivieren.",
		"localStorageNotAvailable": "Ihr Browser scheint zu alt zu sein, um diese Einstellung zu unterstützen. Bitte bringen Sie Ihren Browser auf den neuesten Stand für mehr Kontrolle über Ihre Daten."
	}
};

/**
 * display the current status of the tracking (enabled or not)
 */
function matomoDisplayStatus()
{
	$('.matomo-optout').each(function() {
		$(this).find('.js').css('display', 'inline');
		$(this).find('.nojs').css('display', 'none');

		if (localStorage.getItem('matomoTrackingEnabled') === 'true') {
			$(this).find('input[name="matomo-optout"]').prop('checked', true);
			for (var langAct in langSnippets) {
				$(this).find('label[for="matomo-optout-' + langAct + '"]').html(langSnippets[langAct].trackingActive);
			}
		}
		if (localStorage.getItem('matomoTrackingEnabled') === 'false') {
			$(this).find('input[name="matomo-optout"]').prop('checked', false);
			for (var langInact in langSnippets) {
				$(this).find('label[for="matomo-optout-' + langInact + '"]').html(langSnippets[langInact].trackingInActive);
			}
		}
	});
}

/**
 * change the status of the tracking
 * called on checkbox click
 */
function matomoChangeStatus()
{
	if (typeof(Storage) !== 'undefined') {
		localStorage.matomoTrackingEnabled = (localStorage.getItem('matomoTrackingEnabled') === 'true') ? 'false' : 'true';
		matomoDisplayStatus();
	}
}

/**
 * get the browser’s DoNotTrack setting
 */
var dnt = (navigator.doNotTrack === "yes" || navigator.doNotTrack === "1" || navigator.msDoNotTrack === "1" || window.doNotTrack === "1") ? true : false;

if (dnt && doNotTrackRespected) {
	/**
	 * if browser DoNotTrack setting is activated show doNotTrackActive text
	 */
	for (var langDNT in langSnippets) {
		$('.matomo-optout[lang="' + langDNT + '"]').html(langSnippets[langDNT].doNotTrackActive);
	}
} else {
	if (typeof(Storage) !== 'undefined') {
		/**
		 * if localStorage exists show status text and set event listener to checkbox
		 */
		matomoDisplayStatus();
		$('input[name="matomo-optout"]').on('change', function() {
			matomoChangeStatus();
		});
	} else {
		for (var langOldBrowser in langSnippets) {
			$('.matomo-optout[lang="' + langOldBrowser + '"]').html(langSnippets[langOldBrowser].localStorageNotAvailable);
		}
	}
}

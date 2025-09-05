---
layout: archive
title: "Notes & Blurbs"
permalink: /blurbs/
author_profile: true
collection: blurbs
---

{% include base_path %}

On this page, you will find my scribbles and notes from coursework, research, or personal projects, mostly on topics that interest me or that I found helpful while studying. Since my mathematical training comes from a hodgepodge of courses that I took and audited, as well as other (mostly serious) reading materials, I cannot guarantee that these blurbs are satisfactory to the rigor police despite my best effort. Any feedback, suggestions, and/or corrections are most welcome! 

> Aside from projects that are in progress, I am still working on transcribing many old notes to make them available here, and this statement is true until I remove this notice. 

{% for post in site.blurbs reversed %}
  {% include archive-single.html %}
{% endfor %}
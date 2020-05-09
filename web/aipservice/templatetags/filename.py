import os
from django import template

register = template.Library()

@register.filter
def filename(filepath):
    result = None
    if filepath:
        result = os.path.basename(filepath)
    return result

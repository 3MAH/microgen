{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :no-members:
   :no-inherited-members:
   :no-special-members:

   .. automethod:: __init__

  {% block methods %}  
  {% if methods %}

    **Methods**
    
      .. autosummary::
         :toctree:
         
      {% for item in all_methods %}
         {%- if not item.startswith('_') or item in ['__call__', '__mul__', '__getitem__', '__len__'] %}
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
      {% for item in inherited_members %}
         {%- if item in ['__call__', '__mul__', '__getitem__', '__len__'] %}
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
      
  {% endif %}
  {% endblock %}

  {% block attributes %}
  {% if attributes %}
  
      .. autosummary::
         :toctree:
         
      {% for item in all_attributes %}
         {%- if not item.startswith('_') %}
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
  {% endif %}
  {% endblock %}
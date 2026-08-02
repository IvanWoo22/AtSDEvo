-- Allow XeLaTeX to wrap Chinese text without requiring the xeCJK package.
local function is_cjk(codepoint)
  return
    (codepoint >= 0x2E80 and codepoint <= 0x2FFF) or
    (codepoint >= 0x3000 and codepoint <= 0x303F) or
    (codepoint >= 0x3400 and codepoint <= 0x4DBF) or
    (codepoint >= 0x4E00 and codepoint <= 0x9FFF) or
    (codepoint >= 0xF900 and codepoint <= 0xFAFF) or
    (codepoint >= 0xFF00 and codepoint <= 0xFFEF)
end

function Str(element)
  local output = pandoc.List()
  local latin_buffer = ""

  local function flush_latin()
    if latin_buffer ~= "" then
      output:insert(pandoc.Str(latin_buffer))
      latin_buffer = ""
    end
  end

  for _, codepoint in utf8.codes(element.text) do
    local character = utf8.char(codepoint)
    if is_cjk(codepoint) then
      flush_latin()
      output:insert(pandoc.Str(character))
      output:insert(pandoc.RawInline("latex", "\\allowbreak{}"))
    else
      latin_buffer = latin_buffer .. character
    end
  end
  flush_latin()
  return output
end
